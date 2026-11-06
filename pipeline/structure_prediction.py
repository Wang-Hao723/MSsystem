"""Validated query-side optional model loader; the accepted library source is unchanged."""
import logging
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import threading
import math

torch = None
_LOAD_LOCK = threading.RLock()
_loaded_path = None

logger = logging.getLogger(__name__)

_MODEL_ROOT = Path("/Users/zizi/Desktop/MsQCModel_v2")
_CHECKPOINT_PATH = _MODEL_ROOT / "checkpoints" / "best_model.pt"

_model = None
_vocab = None
_config = None
_device = None


def _ensure_imports():
    msqv2 = str(_MODEL_ROOT)
    if msqv2 not in sys.path:
        sys.path.insert(0, msqv2)


def get_device():
    if torch.backends.mps.is_available():
        return torch.device("mps")
    elif torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


def _load_model_locked(checkpoint_path=None):
    """加载训练好的 MS2SELFIESModel"""
    global _model, _vocab, _config, _device, _loaded_path, torch
    ckpt_path = (Path(checkpoint_path) if checkpoint_path else _CHECKPOINT_PATH).resolve()
    if _model is not None and _loaded_path == ckpt_path:
        return _model, _vocab, _config, _device

    if not ckpt_path.is_file():
        return None, None, None, None
    for required in ('selfies_utils.py', 'model.py', 'utils.py'):
        if not (_MODEL_ROOT / required).is_file():
            raise FileNotFoundError(f"Structure model source not found: {_MODEL_ROOT / required}")
    import torch as torch_runtime
    torch = torch_runtime
    _ensure_imports()
    from selfies_utils import SelfiesVocabulary
    from model import MS2SELFIESModel

    device = get_device()
    logger.info(f"Loading model from {ckpt_path} on {device}")

    ckpt = torch.load(str(ckpt_path), map_location="cpu", weights_only=True)
    config = ckpt.get("config", {})

    # Rebuild vocabulary
    vocab_data = ckpt.get("vocab", {})
    vocab = SelfiesVocabulary()
    if not vocab_data or not vocab_data.get("token_to_idx") or not vocab_data.get("idx_to_token"):
        raise ValueError("Structure checkpoint lacks explicit vocabulary")
    if vocab_data:
        vocab.token_to_idx = {k: int(v) for k, v in vocab_data["token_to_idx"].items()}
        vocab.idx_to_token = {int(k): v for k, v in vocab_data["idx_to_token"].items()}

    # Extract model config from checkpoint (NOT hardcoded)
    model_cfg = config.get("model", {})
    enc_cfg = model_cfg.get("encoder", {})
    dec_cfg = model_cfg.get("decoder", {})
    gnn_cfg = model_cfg.get("gnn", {})

    model = MS2SELFIESModel(
        vocab_size=len(vocab),
        encoder_type=enc_cfg.get("type", "mlp"),
        input_dim=enc_cfg.get("input_dim", 1000),
        latent_dim=enc_cfg.get("latent_dim", 256),
        d_model=dec_cfg.get("d_model", 256),
        nhead=dec_cfg.get("nhead", 4),
        num_decoder_layers=dec_cfg.get("num_layers", 4),
        dim_feedforward=dec_cfg.get("dim_feedforward", 512),
        max_selfies_len=dec_cfg.get("max_len", 128),
        gnn_hidden_dim=gnn_cfg.get("hidden_dim", 128),
        gnn_output_dim=gnn_cfg.get("output_dim", 256),
        gnn_num_layers=gnn_cfg.get("num_layers", 3),
        use_gnn=False,
        use_qc=False,
    ).to(device)

    # Load state dict, ignoring extra keys
    state_dict = ckpt["model_state_dict"]
    # Remove extra keys not in the model
    model_keys = set(model.state_dict().keys())
    filtered_sd = {k: v for k, v in state_dict.items() if k in model_keys}
    skipped = [k for k in state_dict if k not in model_keys]
    if skipped:
        logger.debug(f"Skipped keys: {skipped[:5]}...")

    missing = model_keys - set(filtered_sd)
    if missing:
        raise ValueError("Structure checkpoint missing required weights: " + ", ".join(sorted(missing)[:10]))
    model.load_state_dict(filtered_sd, strict=True)
    model.eval()

    logger.info(f"Model loaded (epoch {ckpt.get('epoch', '?')}, vocab={len(vocab)} tokens)")
    _model, _vocab, _config, _device = model, vocab, config, device
    _loaded_path = ckpt_path
    return model, vocab, config, device


def load_model(checkpoint_path=None):
    with _LOAD_LOCK:
        return _load_model_locked(checkpoint_path)


def predict_from_spectrum(peaks, max_mz=1000.0, bin_size=1.0, temperature=1.0):
    """从 MS2 峰列表预测化合物结构"""
    for name, value in (("max_mz", max_mz), ("bin_size", bin_size), ("temperature", temperature)):
        if isinstance(value, bool) or not math.isfinite(float(value)) or float(value) <= 0:
            raise ValueError(f"{name} must be finite and positive")
    model, vocab, config, device = load_model()
    if model is None:
        return None

    _ensure_imports()
    import importlib.util
    # 直接从 MsQCModel_v2/utils.py 加载 bin_spectrum（避免与本地 utils 冲突）
    msq_utils_path = _MODEL_ROOT / "utils.py"
    spec = importlib.util.spec_from_file_location("msq_utils", msq_utils_path)
    msq_utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(msq_utils)
    bin_spectrum = msq_utils.bin_spectrum
    
    from selfies_utils import selfies_to_smiles, is_valid_selfies

    valid_peaks = []
    for peak in (() if peaks is None else peaks):
        try:
            mz, intensity = float(peak[0]), float(peak[1])
        except (TypeError, ValueError, IndexError):
            continue
        if math.isfinite(mz) and math.isfinite(intensity) and 0 <= mz < max_mz and intensity > 0:
            valid_peaks.append((mz, intensity))
    if not valid_peaks:
        return None
    mz_list = [p[0] for p in valid_peaks]
    int_list = [p[1] for p in valid_peaks]

    peak_vec = bin_spectrum(mz_list, int_list, max_mz=max_mz, bin_size=bin_size)
    expected_bins = int(config.get("model", {}).get("encoder", {}).get("input_dim", 1000))
    if np.asarray(peak_vec).shape != (expected_bins,) or not np.isfinite(peak_vec).all():
        raise ValueError("Structure input bins do not match the checkpoint or contain nonfinite values")
    peak_tensor = torch.tensor(peak_vec, dtype=torch.float).unsqueeze(0).to(device)

    with torch.no_grad():
        gen_tokens = model.generate(
            peak_tensor, None,
            sos_idx=vocab.sos_idx,
            eos_idx=vocab.eos_idx,
            max_len=model.max_selfies_len,
            temperature=temperature,
        )

    pred_selfies = vocab.decode(gen_tokens[0].tolist())
    pred_smiles = selfies_to_smiles(pred_selfies)
    valid = is_valid_selfies(pred_selfies)

    return {
        "selfies": pred_selfies,
        "smiles": pred_smiles,
        "valid": valid,
        "tokens": gen_tokens[0].tolist(),
    }


def structure_similarity(pred_smiles, candidate_smiles):
    """Morgan fingerprint Tanimoto 相似度 (0-1)"""
    if not pred_smiles or not candidate_smiles:
        return 0.0
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import DataStructs

        mol1 = Chem.MolFromSmiles(pred_smiles)
        mol2 = Chem.MolFromSmiles(candidate_smiles)
        if mol1 is None or mol2 is None:
            return 0.0

        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
        return DataStructs.TanimotoSimilarity(fp1, fp2)
    except Exception:
        return 0.0


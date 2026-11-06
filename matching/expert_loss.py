"""
专家引导的损失函数 — Expert-Guided Loss Functions
==================================================

将质谱裂解原理编码为可微损失函数，在模型训练中作为辅助目标或正则化项。

核心思想:
- 不仅要求模型预测与真值接近(传统 MSE/CrossEntropy)
- 还要求预测结果在化学上合理(符合裂解规则)
- 通过多任务学习的范式，让模型内化专家知识

支持的损失:
1. **FragmentRuleLoss** — 预测碎片需符合已知裂解规则
2. **NeutralLossConsistencyLoss** — 碎片间质量差应对应合理中性丢失
3. **BondEnergyRegularizer** — 预测的断键位点应倾向于低键能位置
4. **FunctionalGroupSignatureLoss** — 特定官能团需产生特征碎片模式
5. **StructureFragmentAlignmentLoss** — 碎片应与预测结构兼容

Author: Expert System Module
Version: 2.0.0
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from .fragmentation_rules import (
    COMMON_NEUTRAL_LOSSES,
    BOND_DISSOCIATION_ENERGIES,
    FRAGMENTATION_RULES,
    FUNCTIONAL_GROUP_SIGNATURES,
    get_neutral_loss_by_mass,
    get_fg_signature,
    estimate_bond_lability,
)

logger = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────
# 辅助: 构建专家知识张量
# ─────────────────────────────────────────────────────────────────


def _build_neutral_loss_prior(
    max_mass: float = 500.0, bin_size: float = 1.0
) -> torch.Tensor:
    """
    构建中性丢失先验向量

    对每个可能的中性丢失质量,编码其"化学优先级"。
    已知的、常见的中性丢失 → 低惩罚
    罕见的、不合理的中性丢失 → 高惩罚
    """
    num_bins = int(max_mass / bin_size)
    prior = torch.zeros(num_bins)

    for nl in COMMON_NEUTRAL_LOSSES.values():
        bin_idx = int(nl.mass / bin_size)
        if 0 <= bin_idx < num_bins:
            # 高优先级 = 高得分(低惩罚)
            prior[bin_idx] = max(prior[bin_idx], nl.priority / 10.0)

    return prior


def _build_bond_energy_matrix(
    num_bond_types: int = 20
) -> torch.Tensor:
    """
    构建键能先验矩阵

    用于为注意力机制提供偏置：
    模型在预测碎片时应"注意"低键能位置
    """
    # 将BDE映射为裂解倾向(归一化到 0-1)
    bde_values = list(BOND_DISSOCIATION_ENERGIES.values())
    if not bde_values:
        return torch.eye(num_bond_types)

    # 按能量排序,最低能量 → 最高倾向
    sorted_bdes = sorted(bde_values, key=lambda x: x.energy_kj_per_mol)
    energies = np.array([b.energy_kj_per_mol for b in sorted_bdes])

    # 归一化
    min_e, max_e = energies.min(), energies.max()
    if max_e > min_e:
        lability = 1.0 - (energies - min_e) / (max_e - min_e)
    else:
        lability = np.ones_like(energies)

    return torch.tensor(lability, dtype=torch.float)


# ─────────────────────────────────────────────────────────────────
# 1. 裂解规则符合度损失
# ─────────────────────────────────────────────────────────────────


class FragmentRuleLoss(nn.Module):
    """
    碎片规则符合度损失

    惩罚无法用已知裂解规则解释的碎片预测。
    将预测的碎片质量分布与已知中性丢失模式进行比较。

    使用 KL 散度或余弦距离：
    L_frag = D_KL(p_pred || p_prior) 或 1 - cos(p_pred, p_prior)
    """

    def __init__(
        self,
        max_mass: float = 500.0,
        bin_size: float = 1.0,
        mode: str = "soft_penalty",
    ):
        """
        Args:
            max_mass: 考虑的最大中性丢失质量
            bin_size: 质量分箱大小
            mode: "soft_penalty" (温和惩罚) 或 "hard_penalty" (强硬惩罚)
        """
        super().__init__()
        self.max_mass = max_mass
        self.bin_size = bin_size
        self.mode = mode

        # 构建先验
        prior = _build_neutral_loss_prior(max_mass, bin_size)
        self.register_buffer("neutral_loss_prior", prior)

    def forward(
        self,
        predicted_fragments: torch.Tensor,  # [batch, num_frags] or [batch, mass_bins]
        precursor_mz: torch.Tensor,          # [batch]
        return_details: bool = False,
    ) -> torch.Tensor:
        """
        计算裂解规则损失

        Args:
            predicted_fragments: 预测的碎片质量或谱图
            precursor_mz: 前体离子质量
            return_details: 是否返回详细信息

        Returns:
            loss 值
        """
        batch_size = predicted_fragments.shape[0]

        if predicted_fragments.dim() == 2:
            # If predicted_fragments is [batch, mass_bins], it's a binned spectrum
            num_bins = predicted_fragments.shape[1]
            if num_bins == len(self.neutral_loss_prior):
                prior = self.neutral_loss_prior.unsqueeze(0).expand(batch_size, -1)
            else:
                # Resize prior
                prior = F.interpolate(
                    self.neutral_loss_prior.view(1, 1, -1),
                    size=num_bins, mode="linear", align_corners=False
                ).squeeze()
        else:
            # [batch, num_frags] — directly fragment mz values
            # Convert to histogram
            num_bins = len(self.neutral_loss_prior)
            frags_hist = torch.zeros(batch_size, num_bins,
                                     device=predicted_fragments.device)
            for b in range(batch_size):
                for mz in predicted_fragments[b]:
                    if mz > 0:
                        diff = abs(precursor_mz[b] - mz)
                        bin_idx = min(int(diff / self.bin_size), num_bins - 1)
                        frags_hist[b, bin_idx] += 1
            predicted_fragments = frags_hist
            prior = self.neutral_loss_prior.unsqueeze(0).expand(
                batch_size, -1).to(predicted_fragments.device)

        # 归一化
        pred_norm = F.normalize(predicted_fragments.float(), p=1, dim=1)
        prior_norm = F.normalize(prior.float(), p=1, dim=1)

        if self.mode == "soft_penalty":
            # 余弦距离: 鼓励形状相似
            cos_sim = F.cosine_similarity(pred_norm, prior_norm, dim=1)
            loss = (1.0 - cos_sim).mean()
        elif self.mode == "hard_penalty":
            # KL散度: 更严厉的惩罚
            eps = 1e-8
            loss = (prior_norm * (prior_norm + eps).log()
                    - prior_norm * (pred_norm + eps).log()).sum(dim=1).mean()
        else:
            # 混合惩罚
            cos_sim = F.cosine_similarity(pred_norm, prior_norm, dim=1)
            kl = (prior_norm * (prior_norm + 1e-8).log()
                  - prior_norm * (pred_norm + 1e-8).log()).sum(dim=1)
            loss = ((1.0 - cos_sim) * 0.5 + kl * 0.5).mean()

        if return_details:
            return loss, {"mode": self.mode, "prior_norm": prior_norm}
        return loss


# ─────────────────────────────────────────────────────────────────
# 2. 中性丢失一致性损失
# ─────────────────────────────────────────────────────────────────


class NeutralLossConsistencyLoss(nn.Module):
    """
    中性丢失一致性损失

    确保模型预测的碎片之间的质量差对应于合理的中性丢失。
    对预测结果中不合理的质量差进行惩罚。

    L_nl = Σ w(diff) · penalty(diff)
    其中 w(diff) 是质量差的先验权重
    """

    def __init__(
        self,
        mass_tolerance: float = 0.5,
        penalty_strength: float = 1.0,
    ):
        super().__init__()
        self.mass_tolerance = mass_tolerance
        self.penalty_strength = penalty_strength

        # 为每个常见中性丢失创建权重
        masses = []
        priorities = []
        for nl in COMMON_NEUTRAL_LOSSES.values():
            masses.append(nl.mass)
            priorities.append(nl.priority / 10.0)

        self.register_buffer("nl_masses", torch.tensor(masses, dtype=torch.float))
        self.register_buffer("nl_priorities", torch.tensor(priorities, dtype=torch.float))

    def forward(
        self,
        predicted_fragments: torch.Tensor,  # [batch, num_frags] fragment mz values
        precursor_mz: torch.Tensor,          # [batch]
    ) -> torch.Tensor:
        """
        Args:
            predicted_fragments: 预测碎片 m/z [batch, max_frags]
            precursor_mz: 前体 m/z [batch]

        Returns:
            loss
        """
        batch_size, num_frags = predicted_fragments.shape
        device = predicted_fragments.device

        total_loss = torch.tensor(0.0, device=device)

        for b in range(batch_size):
            pmz = precursor_mz[b]
            frags = predicted_fragments[b][predicted_fragments[b] > 0]

            if len(frags) < 2:
                continue

            # 计算所有碎片与前体的质量差
            diffs = pmz - frags  # [num_frags]

            for diff in diffs:
                if diff <= 0:
                    # 碎片质量不应大于前体 — 硬惩罚
                    total_loss += 0.5 * self.penalty_strength
                    continue

                # 检查是否有合理的中性丢失匹配
                mass_dists = torch.abs(self.nl_masses.to(device) - diff)
                min_dist = mass_dists.min()

                if min_dist > self.mass_tolerance:
                    # 不合理的中性丢失 — 惩罚(有上限)
                    penalty = (min_dist / self.mass_tolerance - 1.0).clamp(min=0, max=10.0) ** 2
                    total_loss += penalty * 0.1 * self.penalty_strength

        return total_loss / max(batch_size, 1)


# ─────────────────────────────────────────────────────────────────
# 3. 键能正则化
# ─────────────────────────────────────────────────────────────────


class BondEnergyRegularizer(nn.Module):
    """
    键能正则化器

    鼓励模型在预测断键位置时倾向于低键能(高裂解倾向)的位置。
    可以作为注意力偏置或输出层的正则化项。

    使用方法:
        1. 作为注意力偏置: attention_bias += bond_energy_bias
        2. 作为损失项: loss += bond_reg(output, attention_weights)
    """

    def __init__(self, reg_strength: float = 0.01):
        super().__init__()
        self.reg_strength = reg_strength

    def forward(
        self,
        attention_weights: torch.Tensor,  # [batch, seq_len, seq_len]
        bond_energy_features: Optional[torch.Tensor] = None,  # [batch, seq_len]
    ) -> torch.Tensor:
        """
        正则化注意力权重,使其偏向低键能位置

        L_bond = reg_strength * Σ(attention · (1 - lability))
        即高注意力在稳定键上时惩罚更大

        Args:
            attention_weights: 模型注意力权重
            bond_energy_features: 键能特征(可选)

        Returns:
            正则化损失
        """
        if bond_energy_features is None:
            # 没有键能信息时,均匀先验
            return torch.tensor(0.0, device=attention_weights.device)

        batch_size = attention_weights.shape[0]
        loss = torch.tensor(0.0, device=attention_weights.device)

        for b in range(batch_size):
            attn = attention_weights[b]  # [seq, seq]
            bde = bond_energy_features[b]  # [seq]

            # 每个位置的"稳定性" = 1 - 倾向性
            stability = 1.0 - bde.clamp(0, 1)  # 高键能 → 高稳定性

            # 对高注意力-高稳定性的组合进行惩罚
            # (模型不应在稳定键上分配过多注意力)
            penalty = (attn.mean(dim=0) * stability).sum() / max(attn.shape[0], 1)
            loss += penalty

        return loss * self.reg_strength / max(batch_size, 1)

    def get_attention_bias(
        self,
        smiles: str,
        seq_length: int,
    ) -> torch.Tensor:
        """
        生成注意力偏置矩阵

        对低键能位置给予正偏置(增加注意力)，
        对高键能位置给予负偏置(减少注意力)。

        Args:
            smiles: 化合物SMILES
            seq_length: 序列长度

        Returns:
            bias 矩阵 [1, 1, seq_len, seq_len] (用于Transformer)
        """
        lability = estimate_bond_lability(smiles)
        bonds = lability.get("labile_bonds", [])

        # 默认偏置: 全为1(无偏置)
        bias = torch.ones(seq_length, seq_length)

        # 对易断裂键的位置增加偏置
        for bond_info in bonds:
            idx1 = bond_info.get("bond_idx", -1)
            score = bond_info.get("lability_score", 0.5)
            if 0 <= idx1 < seq_length:
                # 在注意力矩阵中对这些位置给予正偏置
                bias[idx1, :] *= (1.0 + score)
                bias[:, idx1] *= (1.0 + score)

        # 归一化并转换为 log-space bias
        bias = bias / bias.mean()
        bias = bias.log()

        return bias.unsqueeze(0).unsqueeze(0)  # [1, 1, seq, seq]


# ─────────────────────────────────────────────────────────────────
# 4. 官能团特征损失
# ─────────────────────────────────────────────────────────────────


class FunctionalGroupSignatureLoss(nn.Module):
    """
    官能团特征签名损失

    确保特定官能团的化合物产生其特征碎片模式。
    如果预测黄酮但没有糖丢失和rDA碎片 → 惩罚。
    """

    def __init__(
        self,
        signature_weight: float = 0.1,
        diag_ion_tolerance: float = 0.05,
    ):
        super().__init__()
        self.signature_weight = signature_weight
        self.di_tol = diag_ion_tolerance

    def forward(
        self,
        predicted_spectrum: torch.Tensor,  # [batch, mass_bins]
        compound_classification: List[str],  # [batch]
        mass_axis: torch.Tensor,            # [mass_bins] mass values
    ) -> torch.Tensor:
        """
        Args:
            predicted_spectrum: 预测的谱图(分箱后)
            compound_classification: 每批的分类标签
            mass_axis: 质量轴

        Returns:
            loss
        """
        batch_size = len(compound_classification)
        loss = torch.tensor(0.0, device=predicted_spectrum.device)

        for b in range(batch_size):
            cls = compound_classification[b]
            sig = get_fg_signature(self._map_class(cls))

            if not sig:
                continue

            expected_ions = sig.get("diagnostic_ions", [])
            spectrum = predicted_spectrum[b]  # [mass_bins]

            # 对每个预期诊断离子,检查预测谱图是否含有对应峰
            missed_penalty = 0.0
            for expected_mz in expected_ions:
                # 找到最近的质量bin
                dists = torch.abs(mass_axis - expected_mz)
                min_idx = torch.argmin(dists)
                min_dist = dists[min_idx]

                if min_dist <= self.di_tol:
                    # 检查该位置的预测强度
                    predicted_intensity = spectrum[min_idx]
                    if predicted_intensity < 0.01:
                        # 应该有一个峰但没有 → 惩罚
                        missed_penalty += 1.0

            if expected_ions:
                loss += missed_penalty / len(expected_ions) * self.signature_weight

        return loss / max(batch_size, 1)

    @staticmethod
    def _map_class(classification: str) -> str:
        """映射分类到官能团"""
        c = classification.lower()
        mapping = {
            "flavonoid": "flavonoid", "alkaloid": "alkaloid",
            "terpenoid": "terpenoid", "terpene": "terpenoid",
            "phenolic": "phenolic_acid", "glycoside": "glycoside",
            "lipid": "fatty_acid", "coumarin": "coumarin",
            "quinone": "quinone",
        }
        for key, val in mapping.items():
            if key in c:
                return val
        return "unknown"


# ─────────────────────────────────────────────────────────────────
# 5. 结构-碎片对齐损失
# ─────────────────────────────────────────────────────────────────


class StructureFragmentAlignmentLoss(nn.Module):
    """
    结构-碎片对齐损失

    确保模型预测的化学结构与碎片模式一致。
    惩罚预测的结构不能解释观测碎片的匹配。

    用于多任务模型同时预测结构和碎片时的一致性约束。
    """

    def __init__(self, alignment_weight: float = 0.05):
        super().__init__()
        self.alignment_weight = alignment_weight

    def forward(
        self,
        structure_embedding: torch.Tensor,  # [batch, struct_dim]
        fragment_embedding: torch.Tensor,    # [batch, frag_dim]
        projection_matrix: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        计算结构嵌入和碎片嵌入之间的对齐损失

        L_align = ||structure_emb - project(fragment_emb)||²

        Args:
            structure_embedding: 来自结构编码器的嵌入
            fragment_embedding: 来自碎片编码器的嵌入
            projection_matrix: 可学习的投影矩阵(可选)

        Returns:
            alignment loss
        """
        if projection_matrix is not None:
            fragment_embedding = fragment_embedding @ projection_matrix

        # 确保维度一致
        if structure_embedding.shape[-1] != fragment_embedding.shape[-1]:
            # 使用线性投影
            if not hasattr(self, "proj"):
                self.proj = nn.Linear(
                    fragment_embedding.shape[-1],
                    structure_embedding.shape[-1],
                ).to(fragment_embedding.device)
            fragment_embedding = self.proj(fragment_embedding)

        # MSE loss
        loss = F.mse_loss(structure_embedding, fragment_embedding)
        return loss * self.alignment_weight


# ─────────────────────────────────────────────────────────────────
# 综合专家损失
# ─────────────────────────────────────────────────────────────────


class CompositeExpertLoss(nn.Module):
    """
    综合专家损失 — 组合所有专家引导的损失项

    L_expert = w1 * L_frag_rule + w2 * L_nl_consistency
             + w3 * L_bond_energy + w4 * L_fg_signature
             + w5 * L_structure_align

    在训练中的使用方式:
        expert_loss_fn = CompositeExpertLoss()
        task_loss = F.mse_loss(pred, target)  # 原任务损失
        expert_loss = expert_loss_fn(
            predicted_spectrum=pred_spectrum,
            precursor_mz=precursor_mz,
            predicted_fragments=pred_frags,
            classification=classifications,
            structure_emb=struct_emb,
            fragment_emb=frag_emb,
            mass_axis=mass_axis,
        )
        total_loss = task_loss + expert_weight * expert_loss
    """

    def __init__(
        self,
        weights: Optional[Dict[str, float]] = None,
        max_mass: float = 500.0,
        bin_size: float = 1.0,
    ):
        super().__init__()

        # 默认权重
        self.weights = weights or {
            "frag_rule": 0.25,
            "nl_consistency": 0.25,
            "bond_energy": 0.15,
            "fg_signature": 0.20,
            "structure_align": 0.15,
        }

        self.frag_rule_loss = FragmentRuleLoss(max_mass, bin_size)
        self.nl_consistency_loss = NeutralLossConsistencyLoss()
        self.bond_energy_reg = BondEnergyRegularizer()
        self.fg_signature_loss = FunctionalGroupSignatureLoss()
        self.structure_align_loss = StructureFragmentAlignmentLoss()

    def forward(
        self,
        predicted_spectrum: Optional[torch.Tensor] = None,
        precursor_mz: Optional[torch.Tensor] = None,
        predicted_fragments: Optional[torch.Tensor] = None,
        classification: Optional[List[str]] = None,
        structure_emb: Optional[torch.Tensor] = None,
        fragment_emb: Optional[torch.Tensor] = None,
        attention_weights: Optional[torch.Tensor] = None,
        bond_energy_features: Optional[torch.Tensor] = None,
        mass_axis: Optional[torch.Tensor] = None,
    ) -> Tuple[torch.Tensor, Dict[str, float]]:
        """
        计算综合专家损失

        Returns:
            (total_loss, loss_breakdown_dict)
        """
        losses = {}
        device = (
            predicted_spectrum.device if predicted_spectrum is not None
            else predicted_fragments.device if predicted_fragments is not None
            else torch.device("cpu")
        )

        # 1. 裂解规则损失
        if (predicted_spectrum is not None and precursor_mz is not None
                and self.weights.get("frag_rule", 0) > 0):
            losses["frag_rule"] = self.frag_rule_loss(
                predicted_spectrum, precursor_mz
            )
        else:
            losses["frag_rule"] = torch.tensor(0.0, device=device)

        # 2. 中性丢失一致性
        if (predicted_fragments is not None and precursor_mz is not None
                and self.weights.get("nl_consistency", 0) > 0):
            losses["nl_consistency"] = self.nl_consistency_loss(
                predicted_fragments, precursor_mz
            )
        else:
            losses["nl_consistency"] = torch.tensor(0.0, device=device)

        # 3. 键能正则化
        if (attention_weights is not None
                and self.weights.get("bond_energy", 0) > 0):
            losses["bond_energy"] = self.bond_energy_reg(
                attention_weights, bond_energy_features
            )
        else:
            losses["bond_energy"] = torch.tensor(0.0, device=device)

        # 4. 官能团特征
        if (predicted_spectrum is not None and classification is not None
                and mass_axis is not None
                and self.weights.get("fg_signature", 0) > 0):
            losses["fg_signature"] = self.fg_signature_loss(
                predicted_spectrum, classification, mass_axis
            )
        else:
            losses["fg_signature"] = torch.tensor(0.0, device=device)

        # 5. 结构-碎片对齐
        if (structure_emb is not None and fragment_emb is not None
                and self.weights.get("structure_align", 0) > 0):
            losses["structure_align"] = self.structure_align_loss(
                structure_emb, fragment_emb
            )
        else:
            losses["structure_align"] = torch.tensor(0.0, device=device)

        # 加权求和
        total_loss = torch.tensor(0.0, device=device)
        for key, weight in self.weights.items():
            total_loss += weight * losses.get(key, 0.0)

        loss_dict = {k: round(v.item(), 6) for k, v in losses.items()}
        loss_dict["total"] = round(total_loss.item(), 6)

        return total_loss, loss_dict


# ─────────────────────────────────────────────────────────────────
# 便捷函数
# ─────────────────────────────────────────────────────────────────

def create_expert_training_bias(
    smiles_list: List[str],
    seq_length: int,
) -> torch.Tensor:
    """
    为一组SMILES创建专家引导的训练偏置

    用于在 Transformer 注意力层中引入键能先验。

    Args:
        smiles_list: SMILES列表
        seq_length: 序列长度

    Returns:
        attention_bias [batch, 1, seq_len, seq_len]
    """
    bond_reg = BondEnergyRegularizer()
    biases = []
    for smi in smiles_list:
        bias = bond_reg.get_attention_bias(smi, seq_length)
        biases.append(bias)
    return torch.cat(biases, dim=0)

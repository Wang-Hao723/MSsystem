# MSsystem 投稿补充源码 — 2026-11-02

本压缩包包含 V16 当前工作树中经审查的 855 个源码文件。它是当前有未提交改动的工作树快照，不是 Git 提交 `e03b36cc5ebd287df9e273c24403b644999d2daa` 的干净导出。

## 授权与许可证

用户于 2026-11-02 明确回复：“确认有权，按MIT和CC BY 4.0准备投稿补充文件”。该回复把原 MIT 授权扩展到 `CURRENT_LICENSE_SCOPE.json` 所列 855 文件中的全部自研源码。项目自研源码按根目录 `LICENSE` 的 MIT License 提供。

第三方源码和浏览器包没有被重新许可为项目自研代码。Cytoscape 3.28.1、Plotly.js 2.35.2、SmilesDrawer 2.0.3，以及 MassCube 衍生代码与测试夹具继续遵守各自原 MIT 声明和版权归属；详见 `THIRD_PARTY_NOTICES.md`、文件内原始声明和 `CURRENT_LICENSE_SCOPE.json` 的逐文件标记。

同一回复授权以 CC BY 4.0 准备作者自有数据；第三方数据继续遵守原许可。数据许可随独立数据补充包提供，本源码包不改变数据许可。

本次准备没有执行公开上传、远程仓库发布、DOI 注册或存档编号申请；`CITATION.cff` 不含虚构 DOI。

## 环境与安装

记录运行环境为 macOS arm64、CPython 3.13.5。`requirements.lock` 固定直接依赖版本，其中 SciPy 为 1.16.3。旧 `environment.yml` 指定 Python 3.10，仅作为历史项目文件保留，不作为本快照的安装建议。

建议建立隔离环境；本次许可整理没有验证全新安装：

```bash
python3.13 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.lock
```

部分依赖与操作系统和 CPU 架构有关；锁文件只固定直接依赖，不是完整传递依赖导出。

## 外部输入

本源码包不能单独完成科学复现。生产 SQLite 谱图库、映射表、大型修复登记表、模型权重、原始 LC-MS/MS 输入及冻结评估产物另行提供。详见 `DEPENDENCY_EXCLUSIONS.json` 和论文数据补充包。

至少配置兼容的外部共享数据目录和可写运行目录：

```bash
export MSSYSTEM_SHARED_DATA_ROOT=/absolute/path/to/Mssystem_shared_data
export MSSYSTEM_RUNTIME_ROOT=/absolute/path/to/a/writable/runtime
export MSSYSTEM_LIBRARY_WRITE_ENABLED=0
```

`config/settings.py` 中的环境变量可覆盖工作站路径。`config/runtime_manifest.json` 与确认授权前的审查包保持相同净化字节：本机绝对路径替换为 `${MSSYSTEM_LOCAL_PATH}/<文件名>`，修改位置记录在 `SOURCE_PRIVACY_SCAN.json`。

## 启动与 API 示例

安装依赖并取得兼容外部输入后：

```bash
MSSYSTEM_PYTHON="$PWD/.venv/bin/python" bash run.sh --foreground --no-open --port 8516
curl -fsS http://127.0.0.1:8516/api/health
curl -fsS http://127.0.0.1:8516/openapi.json -o openapi.json
curl -fsS -F 'files=@example.mgf' http://127.0.0.1:8516/api/upload
```

上传端点的 multipart 字段名为 `files`，支持 `.mgf`、`.mzML`、`.mzXML` 和系统可解析的原始报告 `.html`。本次源码许可整理未启动服务、加载数据库、上传样本或执行科学分析。

## 验证

当前 855 个源码载荷已与授权确认前压缩包逐字节核对。许可整理没有改动 Python 源码。重建后已检查 ZIP CRC、路径安全、压缩包成员集合和全部逐文件 SHA-256。Python 字节未变，因此没有重复 AST 检查；此前结果继续保留在 `INSTALLATION_VERIFICATION.json`。


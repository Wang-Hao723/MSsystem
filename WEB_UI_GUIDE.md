# MSsystem Web UI 使用指南

## 一、快速启动

### 1. 安装依赖

```bash
pip install -r requirements.txt
```

> **注意**：PDF 导出功能需要 `weasyprint`，macOS 用户需先安装系统依赖：
> ```bash
> brew install cairo pango gdk-pixbuf libffi
> pip install weasyprint
> ```
> 如果不需要 PDF 功能，可以跳过 weasyprint 安装，其他功能不受影响。

### 2. 启动服务

```bash
# 方式一：直接启动
python web/server.py

# 方式二：使用 uvicorn（推荐，支持热重载）
uvicorn web.server:app --host 0.0.0.0 --port 8502 --reload
```

启动成功后终端显示：
```
INFO:     Uvicorn running on http://0.0.0.0:8502
INFO:     MSsystem Server started. N tasks loaded.
```

### 3. 访问地址

- **Web 界面**：浏览器打开 `http://localhost:8502`
- **API 文档**：`http://localhost:8502/docs`（Swagger UI 交互式文档）
- **备用文档**：`http://localhost:8502/redoc`

### 4. 旧版 Streamlit 界面（兼容保留）

```bash
streamlit run web/app.py
# 访问 http://localhost:8501
```

---

## 二、功能模块说明

### 📤 数据上传模块

**入口**：左侧导航栏 → "数据上传"

| 功能 | 说明 |
|------|------|
| 参考质谱库加载 | 选择 Small（154K 谱图）或 Full（878K 谱图），点击"加载数据库" |
| 文件上传 | 支持拖拽上传或点击选择，格式：`.mgf` / `.mzML` / `.mzXML` |
| 批量上传 | 可同时选择多个文件，每个文件自动创建独立任务 |
| 文件验证 | 自动检查格式和文件大小，超过 500MB 会提示 |
| 上传进度 | 实时显示上传百分比进度条 |
| 分析参数 | Cosine 阈值、Final Score 阈值、m/z 容差（Da 和 ppm） |

**操作流程**：
1. 选择数据库模式 → 点击"加载数据库" → 等待加载完成
2. 拖拽或点击上传质谱数据文件
3. 调整分析参数（可选）
4. 点击"🚀 开始分析"
5. 自动跳转到"任务管理"页面

---

### 📋 任务管理模块

**入口**：左侧导航栏 → "任务管理"

| 功能 | 说明 |
|------|------|
| 任务列表 | 显示所有任务：ID、文件名、上传时间、状态、进度、匹配数 |
| 实时进度 | 分析过程中通过 SSE 实时推送进度更新（无需手动刷新） |
| 搜索 | 按任务 ID 或文件名模糊搜索 |
| 状态筛选 | 下拉筛选：全部 / 待处理 / 运行中 / 已完成 / 失败 |
| 操作按钮 | 开始分析 / 重新分析 / 查看结果 / 删除 |
| 通知 | 分析完成时弹出 Toast 提示 + 浏览器通知 |

**任务状态说明**：
- ⏳ **待处理**：已上传，等待启动分析
- 🔄 **运行中**：正在分析，实时显示进度百分比和速度
- ✅ **已完成**：分析完成，可查看和导出结果
- ❌ **失败**：分析出错，显示具体错误原因，可重新分析

---

### 🔬 结果查看模块

**入口**：左侧导航栏 → "结果查看"

| 功能 | 说明 |
|------|------|
| 任务选择 | 下拉选择已完成的任务 |
| 摘要统计 | 总谱图数、匹配数、匹配率、高置信度数量 |
| 匹配列表 | 按 Final Score 降序排列的表格 |
| 详情弹窗 | 点击行查看：化合物信息、评分、结构式、专家分析 |

**详情弹窗包含**：
- 化合物名称、分子式、分类、SMILES、InChIKey
- 评分卡片：Cosine、Final、Mass、Fragment、Intensity、PPM
- 专家分析 8 维度评分条形图 + 化学推理说明
- 2D 化学结构式图片

---

### 📥 结果导出模块

**入口**：左侧导航栏 → "结果导出"

| 导出格式 | 说明 |
|----------|------|
| 📊 CSV 结果表 | 匹配化合物名称、得分、分类 |
| 📈 CSV 报告表 | 完整评分与元数据 |
| 📗 Excel 工作簿 | 3 个 Sheet：Matched / Unmatched / Summary |
| 🌐 HTML 标准报告 | 交互式表格含结构图 |
| 🔬 HTML 详细报告 | 碎片匹配 + 双谱图对比 + 结构式 + 专家分析 |
| 🖨️ **PDF 报告** | 服务端生成完整 PDF（需安装 weasyprint） |
| 🧬 结构式图片 | PNG / SVG 格式化合物结构图 |

**批量导出**：勾选多个已完成任务 → 点击"📦 批量导出 CSV (ZIP)" → 下载打包文件

**搜索**：支持按任务 ID 或文件名检索已完成任务

---

## 三、API 接口参考

完整 API 文档请访问 `http://localhost:8502/docs`

### 核心 API 列表

| 方法 | 路径 | 说明 |
|------|------|------|
| `GET` | `/api/library/status` | 获取数据库加载状态 |
| `POST` | `/api/library/load?db_mode=small\|full` | 加载/重新加载参考库 |
| `POST` | `/api/upload` | 上传文件（multipart/form-data） |
| `GET` | `/api/tasks` | 任务列表（可按状态筛选） |
| `GET` | `/api/tasks/search?q=keyword` | **搜索任务**（按 ID/文件名） |
| `GET` | `/api/tasks/{id}` | 获取任务详情 |
| `GET` | `/api/tasks/{id}/status` | 获取任务状态（轻量级） |
| `POST` | `/api/tasks/{id}/analyze` | 启动分析 |
| `POST` | `/api/tasks/{id}/retry` | 重新分析 |
| `DELETE` | `/api/tasks/{id}` | 删除任务 |
| `GET` | `/api/tasks/{id}/results` | 获取分页结果（可筛选） |
| `GET` | `/api/tasks/{id}/download/csv` | 下载 CSV |
| `GET` | `/api/tasks/{id}/download/excel` | 下载 Excel |
| `GET` | `/api/tasks/{id}/download/html` | 下载 HTML 报告 |
| `GET` | `/api/tasks/{id}/download/pdf` | **下载 PDF 报告** |
| `GET` | `/api/tasks/{id}/download/structure/{index}` | 下载结构式图片 |
| `POST` | `/api/export/batch` | 批量导出 ZIP |
| `GET` | `/api/tasks/events` | **SSE 实时事件流** |

---

## 四、常见问题

### Q: 页面打开空白？
检查服务是否启动：终端应显示 `Uvicorn running on http://0.0.0.0:8502`。确认浏览器访问 `http://localhost:8502`。

### Q: 上传文件报错？
- 确认文件格式为 `.mgf`、`.mzML` 或 `.mzXML`
- 确认文件非空
- 确认 `output/tmp_uploads/` 目录存在且有写权限

### Q: 分析一直显示"待处理"？
- 确认已加载数据库（上传页面应显示 ✅ 绿色状态）
- 点击"开始分析"按钮启动

### Q: PDF 导出报错？
需要安装 weasyprint 及其系统依赖：
```bash
# macOS
brew install cairo pango gdk-pixbuf libffi
pip install weasyprint

# Ubuntu/Debian
sudo apt-get install libcairo2-dev libpango1.0-dev
pip install weasyprint
```

### Q: 如何清除新手引导弹窗？
引导仅首次访问时显示。如需重新查看，打开浏览器控制台执行：
```javascript
localStorage.removeItem('mssystem_guide_seen')
```
然后刷新页面。

### Q: 如何同时运行 Streamlit 和 FastAPI？
两个服务使用不同端口，可以同时运行：
```bash
# 终端 1：FastAPI（端口 8502）
uvicorn web.server:app --port 8502 --reload

# 终端 2：Streamlit（端口 8501）
streamlit run web/app.py --server.port 8501
```

---

## 五、技术架构

```
浏览器 (http://localhost:8502)
  │
  ├── index.html          ← 主入口
  ├── css/style.css       ← 设计系统
  └── js/
      ├── api.js          ← API 封装层 + SSE
      ├── app.js          ← 主控制器
      ├── upload.js       ← 数据上传模块
      ├── tasks.js        ← 任务管理模块
      ├── results.js      ← 结果查看模块
      └── export.js       ← 结果导出模块
          │
          ▼
FastAPI Backend (port 8502)
  ├── /api/*              ← REST API 端点
  ├── /api/tasks/events   ← SSE 实时推送
  └── matching/           ← 核心匹配引擎
```

- **前端**：纯 HTML/CSS/JS，无需 npm 构建
- **后端**：FastAPI + SQLite
- **实时通信**：Server-Sent Events (SSE)，自动降级为轮询
- **PDF 生成**：服务端 WeasyPrint（HTML → PDF）

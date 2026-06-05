# PDF Parser Pipeline

将文献 PDF 从原始文件转化为结构化中文 Markdown 的完整工作流，包含注册、解析、匹配和翻译四个阶段。

所有脚本位于 `PDF_parser/` 目录，**可从任意工作目录调用，无需 `cd`**。

---

## 概览

```
raw_pdfs/           →  register_pdfs.py           →  register_pdf/
                                                        ├── registry.json
                                                        ├── pdf_<hash>.pdf
                                                        └── grobid_cache.json

register_pdf/       →  match_registry_to_screen.py →  registry_matched.csv  [可选]

[需要代理时]        →  ecs_proxy.sh create         →  SOCKS5 隧道
register_pdf/       →  run_mineru_literature.py     →  minerU_results/
                                                        └── pdf_<hash>/
                                                             ├── full.md
                                                             └── *.json

minerU_results/     →  translate_mineru.py          →  full_zh.md
```

---

## 依赖与环境

### Python 依赖

```bash
pip install requests pandas openai
```

### 外部服务

| 服务 | 用途 | 默认地址 |
|------|------|----------|
| GROBID | 从 PDF 提取元数据（DOI/PMID/标题） | `http://localhost:8070` |
| MinerU API | PDF 解析为结构化 Markdown | `https://mineru.net` |
| DeepSeek API | 翻译 Markdown | `https://api.deepseek.com` |

### 环境变量

```bash
export MINERU_API_TOKEN=<your_mineru_token>   # run_mineru_literature.py 必需
export DEEPSEEK_API_KEY=<your_deepseek_key>   # translate_mineru.py 必需
export ALIYUN_CLI=/path/to/aliyun             # ecs_proxy.sh 可选，不设则自动查找 PATH
```

---

## 脚本详解

### 1. `register_pdfs.py` — 注册原始 PDF

将原始 PDF 按 SHA-256 哈希重命名并复制到注册目录，生成 `registry.json`。原始文件不会被移动或删除。重复运行是幂等的（已注册的文件会跳过）。

**用法**

```bash
# 使用默认路径（literature/00_similar_data_paper_review/selected_papers/）
python /path/to/PDF_parser/register_pdfs.py

# 指定自定义目录（推荐：用绝对路径）
python /path/to/PDF_parser/register_pdfs.py \
    --raw-dir /your/project/raw_pdfs \
    --reg-dir /your/project/register_pdf
```

**参数**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--raw-dir` | `literature/00_.../raw_pdfs` | 原始 PDF 所在目录 |
| `--reg-dir` | `literature/00_.../register_pdf` | 注册目录（输出） |

**输出**

```
register_pdf/
├── registry.json          # {paper_id: {original_filename, sha256, registered_at}}
└── pdf_<hash16>.pdf       # 按 SHA-256 前16位命名的副本
```

---

### 2. `match_registry_to_screen.py` — 匹配注册 PDF 与筛选 CSV（可选）

调用 GROBID 提取每个 PDF 的元数据，然后按 DOI → PMID → 标题相似度的优先级匹配文献筛选 CSV，输出带匹配结果的 CSV。GROBID 结果会缓存在 `grobid_cache.json`，重复运行只处理新文件。此步骤为元数据富化，不影响解析和翻译流程。

**前提**：GROBID 服务必须在 `http://localhost:8070` 运行。

**用法**

```bash
python /path/to/PDF_parser/match_registry_to_screen.py \
    --reg-dir /your/project/register_pdf \
    --csv     /your/project/screened/my_screened.csv \
    --out     /your/project/registry_matched.csv

# 强制重新运行 GROBID（忽略缓存）
python /path/to/PDF_parser/match_registry_to_screen.py \
    --reg-dir ... --csv ... --out ... --no-cache
```

**参数**

| 参数 | 必需 | 说明 |
|------|------|------|
| `--reg-dir` | 是 | 包含 `registry.json` 和重命名 PDF 的目录 |
| `--csv` | 是 | 文献筛选 CSV（需含 `doi`、`pmid`、`title` 列） |
| `--out` | 是 | 输出 CSV 路径 |
| `--no-cache` | 否 | 重新运行 GROBID，不使用已缓存结果 |

**输出列**

| 列名 | 说明 |
|------|------|
| `paper_id` | `pdf_<hash>` |
| `original_filename` | 原始文件名 |
| `grobid_title/doi/pmid` | GROBID 提取的元数据 |
| `match_method` | `doi` / `pmid` / `title_sim=0.XX` / `unmatched` |
| `csv_pmid/title/doi` | 筛选 CSV 中匹配行的字段 |
| `similarity_score` | 筛选 CSV 的相似性评分 |

---

### 3. `ecs_proxy.sh` — 创建/销毁阿里云 ECS SOCKS5 代理

MinerU 上传文件需要经由阿里云 OSS，从境外服务器访问时需要中国大陆代理节点。此脚本在杭州区域按量付费创建一台最小型 ECS 实例（`ecs.t6-c1m1.large`），用作 SSH SOCKS5 隧道，完成后立即销毁以停止计费。

状态文件 `.ecs_instance_id` 保存在脚本所在的 `PDF_parser/` 目录，**可从任意目录调用，`create` 和 `delete` 无需在同一工作目录执行**。

**前提**：已配置 `aliyun` CLI 并完成鉴权（`aliyun configure`）。

**创建代理**

```bash
/path/to/PDF_parser/ecs_proxy.sh create
```

脚本自动完成：查找/创建 VPC、VSwitch、安全组 → 查找 Ubuntu 22.04 镜像 → 创建实例 → 分配公网 IP → 启动实例，然后打印后续操作指令。

**完整工作流**

```bash
# 终端 A — 建立 SSH 隧道（脚本输出密码后执行）
ssh -D 1080 -N -o ServerAliveInterval=60 root@<PUBLIC_IP>

# 终端 B — 验证隧道可用（应返回 "CN"）
curl --socks5 socks5://127.0.0.1:1080 https://ipinfo.io

# 终端 B — 通过代理运行 MinerU
export ALL_PROXY=socks5://127.0.0.1:1080
python /path/to/PDF_parser/run_mineru_literature.py --skip-existing

# 完成后立即销毁实例（停止计费）
/path/to/PDF_parser/ecs_proxy.sh delete
```

**销毁实例**

```bash
/path/to/PDF_parser/ecs_proxy.sh delete
```

停止并强制删除实例，自动清理 `.ecs_instance_id`。**务必在完成解析后执行**，否则持续计费。

---

### 4. `run_mineru_literature.py` — 上传 PDF 至 MinerU 解析

将注册目录中的 PDF 逐一上传至 MinerU 云端 API（`vlm` 模型），轮询解析结果，下载并解压输出的 ZIP，保留 `.md` 和 `.json` 文件。支持断点续传（通过 `.batch_id` 文件恢复）。

**前提**：`MINERU_API_TOKEN` 环境变量已设置。如需翻墙，先建立 SOCKS5 隧道并设置 `ALL_PROXY`。

**用法**

```bash
# 仅处理第一个 PDF（测试连通性）
python /path/to/PDF_parser/run_mineru_literature.py --test

# 处理所有 PDF（跳过已完成的）
python /path/to/PDF_parser/run_mineru_literature.py --skip-existing

# 同时下载图片（默认只下载 .md 和 .json）
python /path/to/PDF_parser/run_mineru_literature.py --skip-existing --image

# 自定义注册目录和输出目录
python /path/to/PDF_parser/run_mineru_literature.py \
    --reg-dir /your/project/register_pdf \
    --out-dir /your/project/minerU_results
```

**参数**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--test` | — | 仅处理第一个 PDF |
| `--skip-existing` | — | 跳过已有输出目录（无 `.batch_id` 文件）的 PDF |
| `--image` | — | 额外下载图片文件（png/jpg 等） |
| `--reg-dir` | `literature/00_.../register_pdf` | 注册目录 |
| `--out-dir` | `literature/00_.../minerU_results` | 结果输出目录 |

**输出结构**

```
minerU_results/
└── pdf_<hash>/
    ├── full.md        # 完整解析后的 Markdown
    ├── *.json         # 结构化内容（段落、表格等）
    └── images/        # 仅 --image 时包含
```

**断点续传**：处理中断时，每个 PDF 目录下的 `.batch_id` 文件保存了云端任务 ID，下次运行会自动恢复，无需重新上传。

---

### 5. `translate_mineru.py` — 翻译 MinerU 输出的 Markdown

将 MinerU 输出的 `full.md` 翻译为中文（或其他语言），使用 DeepSeek API 并发翻译，按顶级 `#` 标题分段处理，保持 Markdown 结构完整。基因符号、药物名、代码块、LaTeX 数学公式等不会被翻译。

**前提**：`DEEPSEEK_API_KEY` 环境变量已设置。

**用法**

```bash
# 翻译单个文件
python /path/to/PDF_parser/translate_mineru.py \
    --input /your/project/minerU_results/pdf_xxx/full.md

# 翻译整个 minerU_results 目录（自动跳过已翻译）
python /path/to/PDF_parser/translate_mineru.py \
    --dir /your/project/minerU_results/

# 强制覆盖已有翻译
python /path/to/PDF_parser/translate_mineru.py --input .../full.md --overwrite

# 翻译为日语
python /path/to/PDF_parser/translate_mineru.py --input .../full.md --lang Japanese

# 调整并发 API 调用数
python /path/to/PDF_parser/translate_mineru.py --dir .../minerU_results/ --workers 8
```

**参数**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--input FILE` | — | 单个 `full.md` 路径（与 `--dir` 二选一） |
| `--dir DIR` | — | `minerU_results` 目录，自动查找所有 `full.md` |
| `--lang` | `Chinese` | 目标语言（Chinese/Japanese/Korean/French/German/Spanish） |
| `--overwrite` | — | 覆盖已存在的翻译文件 |
| `--workers` | `5` | 每个文件的并发 API 调用数 |

**输出**：在 `full.md` 同级目录生成 `full_zh.md`（中文）或 `full_<lang>.md`。

---

## 完整工作流示例

以下以另一个项目目录 `/your/project/` 为例，**无需 `cd` 进入 PDF_parser 所在的 repo**：

```bash
PDF_PARSER=/ShangGaoAIProjects/Gastric/genome/PDF_parser

# 0. 设置环境变量
export MINERU_API_TOKEN=your_token
export DEEPSEEK_API_KEY=your_key

# 1. 注册 PDF
python $PDF_PARSER/register_pdfs.py \
    --raw-dir /your/project/raw_pdfs \
    --reg-dir /your/project/register_pdf

# 2. （可选）匹配筛选 CSV
python $PDF_PARSER/match_registry_to_screen.py \
    --reg-dir /your/project/register_pdf \
    --csv     /your/project/screened/my_screened.csv \
    --out     /your/project/registry_matched.csv

# 3. 如果需要代理访问 MinerU（境外服务器）
$PDF_PARSER/ecs_proxy.sh create
# 在另一个终端建立 SSH 隧道后：
export ALL_PROXY=socks5://127.0.0.1:1080

# 4. 先测试一个 PDF，确认 MinerU 连通
python $PDF_PARSER/run_mineru_literature.py \
    --reg-dir /your/project/register_pdf \
    --out-dir /your/project/minerU_results \
    --test

# 5. 处理全部 PDF
python $PDF_PARSER/run_mineru_literature.py \
    --reg-dir /your/project/register_pdf \
    --out-dir /your/project/minerU_results \
    --skip-existing

# 6. 销毁代理（停止计费）
$PDF_PARSER/ecs_proxy.sh delete

# 7. 翻译所有 full.md 为中文
python $PDF_PARSER/translate_mineru.py \
    --dir /your/project/minerU_results/
```

---

## 常见问题

**GROBID 返回错误或超时**

确认 GROBID 服务正在运行：
```bash
curl http://localhost:8070/api/isalive
```
应返回 `true`。GROBID 通过 Docker 启动：
```bash
docker ps | grep grobid
```

**MinerU 上传卡在 `waiting-file`**

说明 OSS 上传未完成，通常是网络问题。删除对应的 `.batch_id` 文件后重新运行，或启用 SOCKS5 代理再试。

**`aliyun` 命令找不到**

设置环境变量指向实际路径：
```bash
export ALIYUN_CLI=/path/to/your/aliyun
/path/to/PDF_parser/ecs_proxy.sh create
```

**翻译文件乱序**

`translate_mineru.py` 并发翻译各 section，最终按原始顺序写入，不会乱序。如怀疑内容问题，用 `--overwrite` 重新翻译。

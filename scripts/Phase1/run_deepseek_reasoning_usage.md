# run_deepseek_reasoning.py — Usage Guide

Script location: `scripts/Phase1/run_deepseek_reasoning.py`

Sends a Markdown prompt file to the DeepSeek API and saves the model's reasoning chain and final answer as separate Markdown files.

---

## Requirements

| Requirement | Detail |
|-------------|--------|
| Python env | `micromamba run -n ProjectGeneration` |
| Package | `openai` (already installed) |
| API key | `DEEPSEEK_API_KEY` environment variable (already set in shell) |

---

## Models

| Model ID | Description | Use case |
|----------|-------------|----------|
| `deepseek-reasoner` | DeepSeek-R1 — chain-of-thought reasoning **(default)** | Complex biological/analytical questions |
| `deepseek-chat` | DeepSeek-V3 — no thinking chain | Quick summarisation, formatting tasks |

---

## Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--input` | Yes | — | Path to the prompt `.md` file |
| `--output-dir` | No | Same directory as `--input` | Directory where output files are written |
| `--model` | No | `deepseek-reasoner` | DeepSeek model ID |
| `--system` | No | Oncology/genomics expert prompt (see below) | System prompt override |
| `--no-stream` | No | Streaming enabled | Disable streaming; wait for full response before printing |
| `--no-thinking` | No | Thinking chain saved | Skip saving the `_thinking.md` output file |

### Default system prompt

```
You are an expert computational biologist and oncologist specialising in
gastric cancer genomics, tumour molecular subtyping, and somatic mutation
analysis. Provide mechanistic, evidence-based answers grounded in the data
provided. Where relevant, cite known gastric cancer molecular subtypes
(TCGA: EBV, MSI, GS, CIN; ACRG: MSI, MSS/EMT, MSS/TP53+, MSS/TP53−).
Be precise, structured, and scientific.
```

---

## Output Files

Two files are written, named after the input prompt file stem:

| File | Written when | Content |
|------|-------------|---------|
| `<stem>_thinking.md` | `deepseek-reasoner` model + `--no-thinking` not set | Internal reasoning chain |
| `<stem>_response.md` | Always | Final answer |

Both files include a header with model name, prompt filename, and timestamp.

---

## Usage Examples

### 1. Basic — run with all defaults

```bash
/home/sgao30/micromamba/bin/micromamba run -n ProjectGeneration python3 \
    scripts/Phase1/run_deepseek_reasoning.py \
    --input analysis_results/Phase1/TP53_mutual_exclusive/analysis_prompt.md
```

Output written to:
```
analysis_results/Phase1/TP53_mutual_exclusive/analysis_prompt_thinking.md
analysis_results/Phase1/TP53_mutual_exclusive/analysis_prompt_response.md
```

---

### 2. Send output to a separate directory

```bash
/home/sgao30/micromamba/bin/micromamba run -n ProjectGeneration python3 \
    scripts/Phase1/run_deepseek_reasoning.py \
    --input analysis_results/Phase1/neo4j_queries/results_summary.md \
    --output-dir analysis_results/Phase1/llm_responses/
```

---

### 3. Use the faster chat model (no reasoning chain)

```bash
/home/sgao30/micromamba/bin/micromamba run -n ProjectGeneration python3 \
    scripts/Phase1/run_deepseek_reasoning.py \
    --input my_prompt.md \
    --model deepseek-chat
```

---

### 4. Override the system prompt for a non-oncology task

```bash
/home/sgao30/micromamba/bin/micromamba run -n ProjectGeneration python3 \
    scripts/Phase1/run_deepseek_reasoning.py \
    --input my_prompt.md \
    --system "You are an expert biostatistician. Answer concisely and mathematically."
```

---

### 5. Disable streaming (useful when running in background or redirecting stdout)

```bash
/home/sgao30/micromamba/bin/micromamba run -n ProjectGeneration python3 \
    scripts/Phase1/run_deepseek_reasoning.py \
    --input my_prompt.md \
    --no-stream \
    > run.log 2>&1 &
```

---

### 6. Skip saving the thinking chain (save only the final answer)

```bash
/home/sgao30/micromamba/bin/micromamba run -n ProjectGeneration python3 \
    scripts/Phase1/run_deepseek_reasoning.py \
    --input my_prompt.md \
    --no-thinking
```

---

## How to Write a Good Prompt File

The entire content of the `.md` file is sent as the user message. Structure the file with:

1. **Background** — what dataset, cohort, method was used
2. **Data** — paste or reference the key results (tables, statistics, gene lists)
3. **Specific sub-questions** — numbered questions the model should address
4. **Final instruction** — what kind of answer is expected (mechanistic, statistical, clinical, etc.)

See `analysis_results/Phase1/TP53_mutual_exclusive/analysis_prompt.md` as a reference example.

---

## Tips

- **Long prompts are fine.** The script accepts prompts of any length (the TP53 mutual exclusivity prompt was 17,226 characters and ran successfully).
- **Streaming shows progress in real time.** The reasoning chain is printed first, then the answer. For long reasoning tasks this may take several minutes.
- **`deepseek-reasoner` is billed separately** for reasoning tokens and output tokens. Use `deepseek-chat` for simple formatting or summarisation tasks to save cost.
- **The thinking chain is often as valuable as the answer.** It reveals the model's assumptions and logical steps, and is useful for debugging prompt quality.
- **API key is loaded from the environment** — no need to pass it explicitly. If the key is missing the script exits with a clear error message.

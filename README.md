# GeneMind

**GeneMind** is a pipeline that leverages large language models (LLMs) to analyze a given group of genes within a specific research context.  

## Installation

GeneMind requires the **UV package manager**. You can install it following the instructions [here](https://docs.astral.sh/uv/#installation).  

Clone the repository and set up the environment:

```bash
git clone https://github.com/ErwannCnl/genemind
cd genemind 
uv python pin 3.11
uv sync
source .venv/bin/activate
```
## Usage
```bash
uv run main.py "your prompt as free text"
```
### Prompt Guidelines:
Make sure your prompt includes:
- Gene list
- Organism
- How the data was obtained
- Research context
- Desired analysis (currently only GSEA is supported)

### Prompt example
> I identified these five genes to be significantly more mutated than expected by chance in my cohort of human brain cancer patients: CCL14, TLR4, TLR2, IL1B. Do a GSEA.
---
When finished, this will create a new pdf file called `gene_report_from_text.pdf` in the current working directory containing the report.

## Features
- Perform gene set analysis using LLMs
- Currently supported analyses: **GSEA**
- Research context: **Cancer research**
- Supported organisms: **Homo sapiens**

---


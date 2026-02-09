# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MASON (**M**ake **A**nti**S**ense **O**ligos **N**ow) is a Flask-based web tool for designing antisense oligonucleotide (ASO) sequences targeting bacterial genes. It calculates sequence attributes (melting temperature, off-targets) for any gene in any bacterium. The public instance runs at https://www.helmholtz-hiri.de/en/datasets/mason/.

## Running the Server

All commands must be run from the `browser/` directory.

```bash
# Set up environment (first time only)
conda env create -n mason_environment --file mason_server.yml
conda activate mason_environment

# Start the Flask dev server (runs on http://127.0.0.1:5000/)
python run.py
```

## Architecture

### Web Layer (Flask)

- **`browser/run.py`** — Entry point, imports and runs the Flask app
- **`browser/pnag/__init__.py`** — Flask app factory; loads config from `browser/config.json`, sets up Bcrypt, Mail
- **`browser/pnag/routs.py`** — All route definitions (URL handlers)
- **`browser/pnag/forms.py`** — WTForms form classes for input validation

### Three Main Tools (exposed as separate web routes)

1. **MASON** (`/start`, `/startauto`) — Designs ASO sequences for target gene(s). Accepts locus tags, ASO length (7-16), mismatch tolerance (0-4), and optional upstream bases. Runs `mason.sh` in a background thread.
2. **Scrambler** (`/scrambler`) — Takes a single ASO sequence, generates 500 shuffled controls via `esl-shuffle`, and evaluates off-targets. Runs `scrambler.sh` in a background thread.
3. **ASO-Checker** (`/ASO_checker`) — Evaluates user-provided ASO sequence(s) (FASTA format supported) for off-target binding. Also runs `scrambler.sh` (with `-m checker` mode).

### Computation Pipeline (shell scripts + helper scripts)

- **`browser/mason.sh`** — Main MASON pipeline: extracts gene regions from GFF/FASTA using bedtools, designs ASOs via `pnag/make_pnas.py`, runs off-target analysis with seqmap, post-processes results
- **`browser/scrambler.sh`** — Scrambler/Checker pipeline: shuffles sequences, runs seqmap off-target analysis, summarizes results differently based on mode (scrambler vs checker)
- **`browser/start.py`** — Python wrapper that launches bash scripts via `subprocess.run` and writes `done.txt` sentinel files when complete

### Background Processing Pattern

All computations run in Python `threading.Thread`. The web UI polls for completion by checking for a `done.txt` file in `pnag/static/data/<timestamp>/`. Results are stored in timestamped directories under `pnag/static/data/`. Old result directories (>30 days) are auto-deleted by the shell scripts.

### Key Helper Scripts (in `browser/pnag/`)

| Script | Language | Purpose |
|--------|----------|---------|
| `make_pnas.py` | Python | Generates candidate ASO sequences from target gene regions |
| `modify_PNAs.py` | Python | Processes user-provided PNA sequences |
| `scrambler_modify_pnas.py` | Python | Processes shuffled scrambler sequences |
| `checker_modify_pnas.py` | Python | Processes checker input sequences |
| `summarize_offtargets.py` | Python | Summarizes off-target results for MASON |
| `summarize_ots_scrambler.py` | Python | Summarizes off-target results for Scrambler |
| `modify_gff.R` | R | Adjusts GFF entries that exceed chromosome boundaries |
| `melting.R` | R | Calculates melting temperatures |
| `make_final_table_mason.R` | R | Generates final MASON result tables |
| `summarize_ots_checker.R` | R | Summarizes off-target results for Checker |

### External Tool Dependencies

The pipeline relies on command-line bioinformatics tools (installed via conda): **bedtools**, **seqmap**, **bioawk**, **esl-shuffle** (from ViennaRNA/Infernal), **Rscript** (with Bioconductor packages including rmelting, Biostrings).

### Preset Organisms

Four preset genomes are bundled in `pnag/static/data/presets/`: E. coli K-12, S. Typhimurium SL1344, C. difficile 630, F. nucleatum ATCC 23726. Users can also upload custom FASTA + GFF files.

## Important Notes

- The shell scripts hardcode a conda activation path (`/home/jakob/miniconda3/...`). This must be updated for local development.
- `config.json` contains the Flask secret key and email credentials — do not commit changes to this file.
- The app must be run from the `browser/` directory because paths are relative (e.g., `./pnag/static/data/`).
- Templates are in `browser/pnag/templates/`, static assets in `browser/pnag/static/`.

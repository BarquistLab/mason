# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running the Server

```bash
cd browser
conda activate mason_environment
python run.py
# Runs at http://127.0.0.1:5000/
```

## Setup

```bash
git clone git@github.com:BarquistLab/mason.git
cd browser
conda env create -n mason_environment --file mason_server.yml
conda activate mason_environment
```

## Architecture

MASON (Make AntiSense Oligos Now) is a Flask web app for designing antisense oligonucleotides (ASOs/PNAs) targeting bacterial genes. It has three tools that share pipeline infrastructure:

1. **MASON** — Designs ASOs from target gene(s). Route: `/start`. Pipeline: `mason.sh` → `make_pnas.py` → seqmap → `summarize_ots.R` → ML prediction → final table.
2. **Scrambler** — Generates scrambled negative controls for a given ASO. Route: `/scrambler`. Pipeline: `scrambler.sh` (esl-shuffle → `scrambler_modify_pnas.py` → seqmap → `summarize_ots_scrambler.py`).
3. **ASO-Checker** — Validates user-provided ASO sequences. Route: `/ASO_checker`. Reuses `scrambler.sh` with `mode="checker"` flag → `checker_modify_pnas.py` → `summarize_ots_checker.R`.

### Execution Flow

All three tools follow the same pattern:
- Flask route (in `routs.py`) validates form input → spawns `threading.Thread` → calls `start.py:run_pipeline()` → runs shell pipeline via `subprocess.run()` → completion signaled by writing `done.txt`
- Results stored in timestamped directories under `browser/pnag/static/data/YYYY_MM_DD_HH_MM_SS/`

### Key Files (all under `browser/`)

| Layer | Files |
|-------|-------|
| Entry point | `run.py` |
| Flask app | `pnag/__init__.py`, `pnag/routs.py`, `pnag/forms.py` |
| Job runner | `start.py` (generic `run_pipeline()`) |
| Bash pipelines | `mason.sh`, `scrambler.sh`, `common_pipeline.sh` (shared functions) |
| Python analysis | `pnag/make_pnas.py`, `pnag/modify_PNAs.py`, `pnag/scrambler_modify_pnas.py`, `pnag/checker_modify_pnas.py`, `pnag/pna_utils.py` |
| R analysis | `pnag/melting.R`, `pnag/summarize_ots.R`, `pnag/summarize_ots_checker.R`, `pnag/r_utils.R` |
| ML model | `pnag/ML_run.py`, `pnag/static/rf_optimized_model_mason.sav` |
| Templates | `pnag/templates/` (Jinja2 + Bootstrap 4) |

### Shared Utility Modules

- **`pna_utils.py`** — `calculate_purine_stats()`, `calculate_self_complementarity()` (uses cdifflib)
- **`common_pipeline.sh`** — `setup_directories()`, `copy_reference_files()`, `extract_gff_transcripts()`, `run_seqmap_and_process_mismatches()`
- **`r_utils.R`** — `calculate_pna_mw()`, `count_offtargets_by_mismatch()`

### Important Implementation Details

- `mason.sh` uses `coord_offset=0`; `scrambler.sh` uses `coord_offset=-32` in the AWK mismatch pipeline
- `mason.sh` renames uploaded files (spaces → underscores via `GFF_NEW`/`FASTA_NEW`); `scrambler.sh` does not
- `done.txt` sentinel paths differ: mason writes to `/../done.txt`, scrambler/checker write to `/done.txt`
- `summarize_ots.R` indexes results by row name; `summarize_ots_checker.R` indexes by ASO column filter
- `checker_modify_pnas.py` uses cdifflib (C-accelerated drop-in replacement for difflib)

### External Bioinformatics Tools (must be in PATH)

bedtools, seqmap, bioawk, cusp/cai (EMBOSS), RNAfold (ViennaRNA), varna, esl-shuffle (Easel/HMMER)

## Git Conventions

- Branch: `mason_2.0` (active development), `main` (stable)
- Do NOT add Co-Authored-By lines to commits

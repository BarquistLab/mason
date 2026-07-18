# <img src="./browser/pnag/static/mason.png" alt="drawing" width="500"/>



Authors: Jakob Jung, Patrick Pfau, Lars Barquist

Date: 03-05-2022

This is the github page of the MASON webserver. MASON (**M**ake **A**nti**S**ense **O**ligos Now) is a user-friendly web-tool that can be used to design ASO sequences for any gene in any bacterium of interest. MASON can be accessed via https://www.helmholtz-hiri.de/en/datasets/mason/.

In addition to generating ASO sequences able to bind a target gene, it calculates sequence attributes such as the melting temperature of ASO-RNA interaction and predicts possible off-targets within the targeted microbes.

## Tools

### MASON
Designs antisense oligonucleotide (ASO) sequences for one or more bacterial genes.

- **ASO design** — generates antisense oligo sequences targeting the start codon and the Shine–Dalgarno region (if found)
- **Flexible target window** — optionally design ASOs a chosen number of bases upstream and/or downstream of the start codon's first base
- **Sequence attributes** — melting temperature, purine content, self-complementarity, self-binding free energy, and molecular weight for every ASO
- **Off-target prediction** — screens the target organism's transcriptome, and optionally the human transcriptome, the human microbiome, or a user-defined essential-gene set
- **RNA secondary structure** — visualizes the target's translation initiation region (VARNA), coloured by base-pair probability
- **MIC prediction (optional)** — a random-forest model predicts the minimum inhibitory concentration of each ASO and ranks them, with per-ASO SHAP force plots explaining each prediction

### Scrambler
Generates scrambled control sequences for a given ASO.

- Takes a single ASO sequence and generates 500 shuffled controls via `esl-shuffle`
- Evaluates off-target binding for each shuffled sequence
- Useful for comparing on-target ASOs against a null distribution of random sequences

### ASO-Checker
Evaluates user-provided ASO sequences for off-target binding.

- Accepts one or more ASO sequences (FASTA format supported)
- Screens for off-targets against the target organism's transcriptome
- Reports off-target counts by mismatch level and binding location

## Quick start

```bash
git clone git@github.com:BarquistLab/mason.git
cd mason/browser
conda env create -n browser --file mason_server.yml
conda activate browser
python run.py
```

Then open http://127.0.0.1:5000/ in your browser.

## Requirements

- Linux with bash, [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) (miniconda recommended), and [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
- All Python, R, and CLI dependencies are installed automatically via `mason_server.yml`
- **VARNA** for RNA structure plots — download the jar from https://varna.lisn.upsaclay.fr/ and create a wrapper script at `/home/jakob/bin/varna`:
  ```bash
  mkdir -p ~/bin
  wget https://varna.lisn.upsaclay.fr/bin/VARNAv3-93.jar -O ~/bin/VARNAv3-93.jar
  echo '#!/bin/bash
  java -cp ~/bin/VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd "$@"' > ~/bin/varna
  chmod +x ~/bin/varna
  ```
- **xvfb** for headless rendering of VARNA plots:
  ```bash
  sudo apt-get install -y xvfb
  ```

See [browser/README.md](./browser/README.md) for detailed installation and setup instructions. The code is published under the MIT license.

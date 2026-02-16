# <img src="./browser/pnag/static/mason.png" alt="drawing" width="500"/>



Authors: Jakob Jung, Patrick Pfau, Lars Barquist

Date: 03-05-2022

This is the github page of the MASON webserver. MASON (**M**ake **A**nti**S**ense **O**ligos Now) is a user-friendly web-tool that can be used to design ASO sequences for any gene in any bacterium of interest. MASON can be accessed via https://www.helmholtz-hiri.de/en/datasets/mason/.

In addition to generating ASO sequences able to bind a target gene, it calculates sequence attributes such as the melting temperature of ASO-RNA interaction and predicts possible off-targets within the targeted microbes.

## Quick start

```bash
git clone git@github.com:BarquistLab/mason.git
cd mason/browser
conda env create --file mason_server.yml
conda activate mason_environment
python run.py
```

Then open http://127.0.0.1:5000/ in your browser.

## Requirements

- Linux with bash, [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html), and [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
- All Python, R, and CLI dependencies are installed automatically via `mason_server.yml`
- **VARNA** (optional) for RNA structure plots — download from https://varna.lisn.upsaclay.fr/ and place on your `PATH`

See [browser/README.md](./browser/README.md) for detailed installation and setup instructions. The code is published under the MIT license.

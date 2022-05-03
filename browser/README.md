# MASON Webserver
<img src="./pnag/static/mason.png" alt="drawing" width="500"/>



Date: 13-04-2020

Author: Jakob Jung, Patrick Pfau

Supervision: Lars Barquist

Software needed: Bash shell, Python (v3.7+) [Anaconda](https://docs.anaconda.com/anaconda/install/linux/) / [Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html), git

This is a basic Webserver for using the MASON algorithm. MASON can be accessed via https://www.helmholtz-hiri.de/en/datasets/mason/. Find below the instructions on how to install the server on a local computer or cluster using a bash shell. 



## Installation on local machine

To install MASON, you need a linux-based bash shell and an installed version of [Anaconda](https://docs.anaconda.com/anaconda/install/linux/) or [Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) and [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git). Start by cloning the git repository by:

```bash
git clone git@github.com:BarquistLab/mason.git
```

 Now, you can enter the directory and install all required packages with conda. For this we can use the mason_environment.yaml file in which packages are stored. This can take some minutes:

```bash
cd browser
conda env create -n mason_environment --file mason_server.yml
conda activate mason_environment
```

That's it. Now you're ready to go to run mason from the command line.



## Run the MASON web-server on a local machine

This Webserver runs with flask. To start it, just type:

```bash
python run.py
```

and an output llooking like this should appear:

```
Debug mode: off
Running on http://127.0.0.1:5000/ (Press CTRL+C to quit)
```

Now open the url http://127.0.0.1:5000/ in your browser and use MASON. 

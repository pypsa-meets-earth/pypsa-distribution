# PyPSA-Distribution

## Development Status: **Under development, contributors are welcome to join**

[![Status Linux](https://github.com/pypsa-meets-earth/pypsa-distribution/actions/workflows/ci-linux.yaml/badge.svg?branch=main&event=push)](https://github.com/pypsa-meets-earth/pypsa-distribution/actions/workflows/ci-linux.yaml)
![Size](https://img.shields.io/github/repo-size/pypsa-meets-earth/pypsa-distribution)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/pypsa-meets-earth/pypsa-distribution/main.svg)](https://results.pre-commit.ci/latest/github/pypsa-meets-earth/pypsa-distribution/main)
[![Discord](https://img.shields.io/discord/911692131440148490?logo=discord)](https://discord.gg/AnuJBk23FU)
[![Google Drive](https://img.shields.io/badge/Google%20Drive-4285F4?style=flat&logo=googledrive&logoColor=white)](https://drive.google.com/drive/folders/1U7fgktbxlaGzWxT2C0-Xv-_ffWCxAKZz)
[![Documentation Status](https://readthedocs.org/projects/pypsa-distribution/badge/?version=latest)](https://pypsa-distribution.readthedocs.io/en/latest/?badge=latest)<!-- [![Status Mac](https://github.com/pypsa-meets-earth/pypsa-earth/actions/workflows/ci-mac.yaml/badge.svg?branch=main&event=push)](https://github.com/pypsa-meets-earth/pypsa-earth/actions/workflows/ci-mac.yaml) -->
<!-- [![Status Windows](https://github.com/pypsa-meets-earth/pypsa-earth/actions/workflows/ci-windows.yaml/badge.svg?branch=main&event=push)](https://github.com/pypsa-meets-earth/pypsa-earth/actions/workflows/ci-windows.yaml) -->
   

**PyPSA-Distribution is a multi-energy model for small scale applications in high spatial and temporal resolution.**
It leverages on the open-source tool PyPSA-Earth and aims at addressing the optimization of small-scale applications.
Currently it focuses on electric off-grid applications as it is at an infant state.


## Get involved

There are multiple ways to get involved and learn more about our work.
Please, checkout the [PyPSA-Earth list](https://github.com/pypsa-meets-earth/pypsa-earth)

The recurrent meeting on PyPSA-Distribution is every second Tuesday at 15:00 AM (CEST time) on Discord, please reach out and join! A calendar invitation may be found [here](https://drive.google.com/file/d/1qGsUOKfoBt3FtXnEXiiLC-H16DiQvyVy/view?usp=drive_link).

## Installation

0. **For Windows users**. Please, that windows is supported using WSL (Windows Subsystem Linux), due [Issue 392 Snakemake](https://github.com/snakemake/snakemake/issues/392)
   Windows users may install the WSL as discussed in [this link](https://code.visualstudio.com/docs/remote/wsl), also including the setup of Visual Studio Code.
   
1. Open your terminal at a location where you want to install pypsa-earth. Type the following in your terminal to download the package from GitHub:

   ```bash
      .../some/path/without/spaces % git clone --recurse-submodules https://github.com/pypsa-meets-earth/pypsa-distribution.git
   ```
2. The python package requirements are curated in the `pypsa-earth/envs/environment.yaml` file.
   The environment can be installed using:

   ```bash
      .../pypsa-distribution % conda env create -f pypsa-earth/envs/environment.yaml
   ```

   If the above takes longer than 30min, you might want to try mamba for faster installation:

   ```bash
     (base) conda install -c conda-forge mamba

     .../pypsa-distribution % mamba env create -f pypsa-earth/envs/environment.yaml
   ```

3. For running the optimization one has to install a solver. Check the [installation](https://pypsa-earth.readthedocs.io/en/latest/installation.html) guideline for more details.

4. To use jupyter lab (new jupyter notebooks) **continue** with the [ipython kernel installation](http://echrislynch.com/2019/02/01/adding-an-environment-to-jupyter-notebooks/) and test if your jupyter lab works:

   ```bash
      .../pypsa-earth % ipython kernel install --user --name=pypsa-earth
      .../pypsa-earth % jupyter lab
   ```
5. Verify or install a java redistribution from the [official website](https://www.oracle.com/java/technologies/downloads/) or equivalent.
   To verify the successfull installation the following code can be tested from bash:

   ```bash
      .../pypsa-distribution % java -version
   ```

   The expected output should resemble the following:

   ```bash
      java version "1.8.0_341"
      Java(TM) SE Runtime Environment (build 1.8.0_341-b10)
      Java HotSpot(TM) 64-Bit Server VM (build 25.341-b10, mixed mode)
   ```
6. Extra python package requirements are satisfied by manually installing them. Currently the ramp package is added to existing pypsa-earth environment by using:

   ```bash
      .../pypsa-distribution % conda activate pypsa-earth
      (pypsa-earth) pip install rampdemand
   ```
   
## Test run

- In the folder open a terminal/command window to be located at this path `~/pypsa-distribution/`
- Activate the environment `conda activate pypsa-earth`
- Run a dryrun of the Snakemake workflow by typing simply in the terminal:
  ```bash
  snakemake -j 1 solve_network -n
  ```

  Remove the -n to do a real run.

## Training

- We recently updated some [hackathon material](https://github.com/pypsa-meets-earth/documentation) for PyPSA-Earth, which are relevant also for PyPSA-Distribution. The hackathon contains jupyter notebooks with exercises. After going through the 1 day theoretical and practical material you should have a suitable coding setup and feel confident about contributing.
- The get a general feeling about the PyPSA functionality, we further recommend going through the [PyPSA](https://github.com/PyPSA/PyPSA/tree/master/examples) and [Atlite](https://github.com/PyPSA/atlite/tree/master/examples) examples.

## Questions and Issues

- We are happy to answer questions and help with issues **if they are public**. Through being public the wider community can benefit from the raised points. Some tips. **Bugs** and **feature requests** should be raised in the [**GitHub Issues**](https://github.com/pypsa-meets-earth/pypsa-distribution/issues/new/choose). **General workflow** or **user questions** as well as discussion points should be posted at the [**GitHub Discussions**](https://github.com/pypsa-meets-earth/pypsa-distribution/discussions/categories/q-a) tab. Happy coding.

<!-- ## Documentation

The documentation is available here: [documentation](https://pypsa-earth.readthedocs.io/en/latest/index.html). -->

## Collaborators


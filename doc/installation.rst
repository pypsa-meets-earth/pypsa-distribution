..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Earth authors

  SPDX-License-Identifier: CC-BY-4.0

.. _installation:

##########################################
Installation
##########################################

The subsequently described installation steps are demonstrated as shell commands, where the path before the ``%`` sign denotes the directory in which the commands following the ``%`` should be entered.



Install/Check Dependencies
===============================

The major software dependencies of PyPSA-Distribution coincides with those of PyPSA-Earth,
described at `this link <https://pypsa-earth.readthedocs.io/en/latest/installation.html#install-dependencies>`_.
In particular, after having verified to have installed ``git``, a ``conda`` distribution (e.g. miniconda), and
an ``Integrated Development Environment (IDE)``, such as Visual Studio Code, you can proceed with the following.

Clone the Repository
====================
.. note::

  **For Windows users**. Please, that windows is supported using WSL (Windows Subsystem Linux), due `Issue 392 Snakemake <https://github.com/snakemake/snakemake/issues/392>`_
  Windows users may install the WSL as discussed in `this link <https://code.visualstudio.com/docs/remote/wsl>`_, also including the setup of Visual Studio Code.
  More information in the short_tutorial section

First of all, clone the `PyPSA-Distribution repository <https://github.com/pypsa-meets-earth/pypsa-distribution/>`_ using the version control system ``git``.
The path to the directory into which the ``git repository`` is cloned, must **not** have any spaces.
If you do not have ``git`` installed, follow installation instructions `here <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_.

.. code:: bash

    /some/other/path % cd /some/path/without/spaces

    /some/path/without/spaces % git clone --recurse-submodules https://github.com/pypsa-meets-earth/pypsa-distribution.git

.. note::

  PyPSA-Distribution depends recoursively on PyPSA-Earth, therefore, option `--recurse-submodules` is needed:
  it automatically downloads the PyPSA-Earth repository as well.


Install the Python Environment
==============================

The python package requirements are curated in the `pypsa-earth/envs/environment.yaml` file.
The environment can be installed using:

.. code:: bash

    .../pypsa-distribution % conda env create -f pypsa-earth/envs/environment.yaml

If the above takes longer than 30min, you might want to try mamba for faster installation:

.. code:: bash

    (base) conda install -c conda-forge mamba

    .../pypsa-distribution % mamba env create -f pypsa-earth/envs/environment.yaml

Extra python package requirements are satisfied by manually installing them. Currently the ramp package is added to existing pypsa-earth environment by using:

.. code:: bash

    .../pypsa-distribution % conda activate pypsa-earth

    (pypsa-earth) pip install rampdemand
   

Install the pre-commit (for developers)
=======================================

We are using the `pre-commit <https://pre-commit.com/>`_ framework to run some checks on the code before committing.
The tool provide automated ways to keep the code clean and consistent.
To make sure to install the tool, please run:

.. code:: bash

    .../pypsa-distribution % conda activate pypsa-earth

    .../pypsa-distribution % pre-commit install


Set Configuration File
================================

PyPSA-Earth has several configuration options that must be specified in a ``config.yaml`` file located in the project directory. An example configuration ``config.distribution.yaml`` is maintained in the repository. More details on the configuration options are in :ref:`config` section.

Note that the ``config.distribution.yaml`` is a difference file with respect to the ``config.pypsa-earth.yaml`` one.

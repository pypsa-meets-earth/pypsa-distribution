..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Earth authors

  SPDX-License-Identifier: CC-BY-4.0

.. _introduction:

##########################################
Introduction
##########################################

PyPSA-Distribution aims at providing a global, open, free, fully transparent, and reproducible energy system model for distribution applications.
It is based on the `PyPSA <https://pypsa.org/>`_ framework and uses the `OpenStreetMap <https://www.openstreetmap.org/>`_ data to create a detailed representation of the distribution energy system.

Examples of possible applications of the model are:

- **Microgrids planning**: The model can be used to plan the expansion of a network of microgrids, e.g. to integrate new renewable energy sources or to electrify heating and transport.
- **Distribution system planning**: The model can be used to plan the expansion of the distribution network, e.g. to integrate new renewable energy sources or to electrify heating and transport.
- **Distribution system operation**: The model can be used to operate the distribution network, e.g. to optimize the dispatch of distributed energy resources or to avoid congestions.
- **Distribution network tariff design**: The model can be used to design tariffs for the distribution network, e.g. to incentivize the use of renewable energy sources or to avoid congestions.
- **Are we missing something?**: Please let us know if you have any other ideas for applications of the model!



Workflow
========

As `PyPSA-Earth <https://pypsa-earth.readthedocs.io/en/latest/>`_,
the generation of the model is controlled by the workflow management system `Snakemake <https://snakemake.bitbucket.io/>`_. In a nutshell,
the ``Snakefile`` declares for each python script in the ``scripts`` directory a rule which describes which files the scripts consume and
produce (their corresponding input and output files). The ``snakemake`` tool then runs the scripts in the correct order according to the
rules' input/output dependencies. Moreover, it is able to track, what parts of the workflow have to be regenerated, when a data file or a
script is modified/updated. Please, refer to PyPSA-Earth documentation for more.


Folder structure
================

The content in this package is organized in folders as described below; for more details, please see the documentation.

- ``data``: Includes input data that is not produced by any ``snakemake`` rule.
- ``scripts``: Includes all the Python scripts executed by the ``snakemake`` rules.
- ``resources``: Stores intermediate results of the PyPSA-Distribution workflow which can be picked up again by subsequent rules.
- ``networks``: Stores intermediate, unsolved stages of the PyPSA network that describes the energy system model of PyPSA-Distribution.
- ``results``: Stores the solved PyPSA network data, summary files and plots.
- ``benchmarks``: Stores ``snakemake`` benchmarks.
- ``logs``: Stores log files about solving, including the solver output, console output and the output of a memory logger.
- ``envs``: Stores the conda environment files to successfully run the workflow.


License
=======

PyPSA-Distribution work is released under multiple licenses, equivalent to those of PyPSA-Earth.

* All original source code is licensed as free software under `GPL-3.0 License <https://github.com/pypsa-meets-earth/pypsa-earth/blob/main/LICENSE>`_.
* The documentation is licensed under `CC-BY-4.0 <https://creativecommons.org/licenses/by/4.0/>`_.
* Configuration files are mostly licensed under `CC0-1.0 <https://creativecommons.org/publicdomain/zero/1.0/>`_.
* Data files are licensed under different licenses as noted below.

For licenses and urls of the data used in PyPSA-Earth, please, refer to `PyPSA-Earth License section <https://pypsa-earth.readthedocs.io/en/latest/introduction.html#license>`_.

* *BY: Attribute Source*
* *NC: Non-Commercial Use Only*
* *SA: Share Alike*

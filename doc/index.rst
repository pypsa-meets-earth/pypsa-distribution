.. PyPSA meets Earth documentation master file, created by
   sphinx-quickstart on Sat May 15 22:52:54 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the PyPSA-Distribution documentation!
================================================

.. note::
    This documentation is under construction and will be updated soon!
    Anyone interested in the project is welcome to join and collaborate with us!

.. image:: https://img.shields.io/github/v/release/pypsa-meets-earth/pypsa-distribution?include_prereleases
    :alt: GitHub release (latest by date including pre-releases)

.. image:: https://github.com/pypsa-meets-earth/pypsa-distribution/actions/workflows/ci-linux.yaml/badge.svg
    :target: https://github.com/pypsa-meets-earth/pypsa-distribution/actions

.. image:: https://readthedocs.org/projects/pypsa-distribution/badge/?version=latest
    :target: https://pypsa-distribution.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://img.shields.io/github/repo-size/pypsa-meets-earth/pypsa-distribution
    :alt: GitHub repo size

.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
    :target: https://www.gnu.org/licenses/gpl-3.0

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
    :alt: Code style Black

.. image:: https://results.pre-commit.ci/badge/github/pypsa-meets-earth/pypsa-distribution/main.svg
    :target: https://results.pre-commit.ci/latest/github/pypsa-meets-earth/pypsa-distribution/main
    :alt: Pre-commit CI-status

.. image:: https://img.shields.io/discord/911692131440148490?logo=discord
    :target: https://discord.gg/AnuJBk23FU
    :alt: Discord

.. image:: https://img.shields.io/badge/Google%20Drive-4285F4?style=flat&logo=googledrive&logoColor=white
    :target: https://drive.google.com/drive/folders/1U7fgktbxlaGzWxT2C0-Xv-_ffWCxAKZz
    :alt: Google Drive

*Motivation*. Closed-source models are the current standard for most policy and industry decisions. However, open models have proven to be
competitive alternatives that promote science, robust technical analysis, collaboration and transparent policy decision making.
Yet, two issues slow the adoption: open models are often designed with limited geographic scope, hindering synergies to collaborate,
or are based on low spatially resolved data, limiting their utility.

*PyPSA-Distribution* aims to promote open-source global energy system model at distribution scale with data in high spatial and temporal resolution.
It enables large-scale collaboration by providing a tool that can model the distribution system of any region in the world.
This work leverages on significant previous work by PyPSA-Earth, PyPSA-Eur and the references and contributions discussed in the repository.
from the European PyPSA-Eur model using new data and functions. It is suitable for operational as well as combined generation,
storage and transmission expansion studies. We work hard to extend the PyPSA-Earth model by end of this year to include sector-coupling,
myopic and perfect pathway expansion capabilities.

Example of desired studies are: microgrids planning, distribution system planning, distribution system operation, distribution network tariff design, ...
**Are we missing something?** Please let us know if you have any other ideas for applications of the model!

*PyPSA meets Earth initiative* members are maintaining the *PyPSA-Distribution* repository.
The `website <https://pypsa-meets-earth.github.io/>`_ provides more context of the initiative and the associated projects. 

==============
Get Involved
==============

Discussions on the PyPSA-Distribution tool are hosted on the `PyPSA meets Earth Discord <https://discord.gg/AnuJBk23FU>`_.

The recurrent meeting on PyPSA-Distribution is every second Tuesday at 17:00 AM (CEST time) on Discord, please reach out and join! A calendar invitation may be loaded `here <https://drive.google.com/file/d/1qGsUOKfoBt3FtXnEXiiLC-H16DiQvyVy/view?usp=drive_link>`_.

=============
Documentation
=============

**Getting Started**

* :doc:`introduction`
* :doc:`installation`
* :doc:`short_tutorial`


.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Getting Started

   introduction
   installation
   short_tutorial

**Model Costumization**

* :doc:`customization_basic`
* :doc:`Run_your_personal_case`
  

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Model Costumization

    customization_basic
    Run_your_personal_case



.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Model Costumization

    customization_basic


**Work flow and API**

* :doc:`api_reference`

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Work flow and API

   api_reference

**Help and References**

* :doc:`release_notes`
* :doc:`how_to_contribute`

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Project Info

   release_notes
   how_to_contribute
   

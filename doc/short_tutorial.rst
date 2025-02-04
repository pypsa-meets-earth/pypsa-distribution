.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _short_tutorial:


##########################################
Short tutorial
##########################################

.. note::

    For realistic case studies it will be necessary to run the tool in full configuration and not in tutorial configuration. The tutorial mode is used only to become familiar with the tool.

The installation procedure installs PyPSA-Distribution model with all the software dependencies needed to build and run it.
The purpose of this section is to make available to the user a simplified case study through which to become familiar with the code and its execution.
To this aim, we have chosen to present the simulation of a microgrid in Nigeria, choosing the tutorial mode within the configuration file. In this way, the user can explore the different possibilities of the tool,
but at a lower expense of computational resources than in the case of a full simulation.

Before we get into the heart of the tutorial, it is important to make a few preliminaries for using the tool if you have a Windows operating system.

Preliminary operations
---------------------------

For windows users: due to a problem with Snakemake, using windows is supported as long as you use WSL, you can easily install it in your device by following these steps.
Open the command prompt and run the following command

.. code:: bash

    wsl --install

Follow the procedure shown in the prompt, and when finished, restart the device.

You may probably need to install conda and create the pypsa-earth environment in the WSL execution environment as well.
To do this you need to first download the linux version from https://docs.anaconda.com/free/miniconda/. 
Once the download is done open the wsl terminal ( you can just type wsl in the search panel), place yourself in the folder where the file you just downloaded is located, 
and run the following commands.

.. code:: bash

    /your-path (base) % bash <conda-installer-name>-latest-Linux-x86_64.sh

Where conda-installer-name will be one of "Miniconda3", "Anaconda", or "Miniforge3".
Follow the prompts on the installer screens. If you are unsure about any setting, accept the defaults. You can change them later.

To check if the installation was successful you can run the following command:

.. code:: bash

    /your-path (base) % conda list

If you see a list of packages after execution, everything should be okay.

To create the environment, simply open the wsl terminal, open to the pypsa-earth folder and run the following commands:

.. code:: bash

    /pypsa-earth % conda activate base

After that:

.. code:: bash

    /pypsa-earth % conda install -c conda-forge mamba

Finally:

.. code:: bash

    /pypsa-earth % mamba env create -f envs/environment.yaml

Now everything should be ready to run the tutorial!

Set the IDE:
---------------------
.. note::

    In this example we will refer to the use of the Visual Studio Code IDE.

First, we have to set up the IDE in order to run the tutorial. 
In particular it will be necessary, first, to open a new terminal. To do this, simply select the "terminal" option 
at the top left of the VSC interface and select new terminal.

Once this is done we will see a window appear at the bottom of the interface. 
For Windows users by default, the terminal will be a cmd, however for proper execution of the tool it is necessary to move to a WSL terminal. 
To do this, it is sufficient to select the down arrow next to "cmd" and select from the drop-down menu the terminal "wsl"

At this point we can move on to activate the previously created pypsa-earth environment ( if you have not already created the pypsa-earth environment you can see how to do it in the "installation" section).
To activate the environment simply run the following command in the terminal you just opened:

.. code:: bash

    .../your-folder (base) % conda activate pypsa-earth

Now you just have to move to the folder in which you want to work, and you'll be ready to get into the thick of the simulation.

.. code:: bash

    .../your-folder (pypsa-earth) % cd your-work-folder


Run the tutorial model
---------------------

A tutorial data kit was developed to facilitate exploring the model.
The user can explore the majority of the model's functions on a local machine by running the tutorial, which uses fewer computational resources than the entire model does. 
Currently, the tutorial case study refers to a microgrid in Nigeria whose coverage area is defined by a rectangle whose vertices have the following coordinates:

-	lon_max: 5.0998
-	lon_min: 6.1700
-	lat_min: 8.2356
-	lat_max: 9.8012

Before actually running the tool, it is always a good idea to check how it will look by using -dryrun or -n Snakemake option:

.. code:: bash

    .../pypsa-distribution (pypsa-earth) % snakemake -j 1 solve_all_networks --dryrun


To run the whole modeling workflow you just need the following command:

.. code:: bash

    .../pypsa-distribution (pypsa-earth) % snakemake -j 1 solve_network

.. note::

    Before running these commands always make sure:
    - you are in the correct folder ( i.e., the folder related to the project where the snakefile is located)
    - that you have enabled the pypsa-earth environment

.. TODO Explain settings of the tutorial case

This command will trigger loading of the whole dataset needed to build the model for a tutorial case if both tutorial and retrieve_databundle flags are on. 
The tutorial model run simulation will take a while (about 20..50 minutes).
If the simulation was successful, you should be able to display an elec.nc file representing 
the optimised network in the networks/results folder.

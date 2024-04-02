.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _short_tutorial:


##########################################
Short tutorial
##########################################

The installation procedure installs PyPSA-Distribution model with all the software dependencies needed to build and run it.
It is first crucial to get familiar with a tutorial where a simpler model is considered. This section explains how to:
-	Run  the tutorial model
-	Run the tool for different case studies than the one in the tutorial

Before we get into the heart of the tutorial, it is important to make a few preliminaries for using the tool if you have a Windows operating system.

Preliminary operations
---------------------------

For windows users: due to a problem with Snakemake using windows is supported as long as you use WSL, you can easily install it in your device by following these steps
Open the command prompt and run the following command

.. code:: bash

    wsl --install

Follow the procedure shown in the prompt, and when finished, restart the device.

You may probably need to install conda and create the pypsa-earth environment in the WSL execution environment as well.
To do this you need to first download the linux version from https://docs.anaconda.com/free/miniconda/. 
Once the download is done open the wsl terminal ( you can just type wsl in the search panel),place yourself in the folder where the file you just downloaded is located,  and run the following commands.

.. code:: bash
    \your-path (base)  $ bash <conda-installer-name>-latest-Linux-x86_64.sh

Where conda-installer-name will be one of "Miniconda3", "Anaconda", or "Miniforge3".
Follow the prompts on the installer screens. If you are unsure about any setting, accept the defaults. You can change them later.

To check if the installation was successful you can run the following command:

.. code:: bash
    \your-path (base)  $ conda list

If you see a list of packages after execution, everything should be okay.

To create the environment, simply open the wsl terminal, open to the pypsa-earth folder and run the following commands:

.. code:: bash
    \pypsa-earth   $ conda activate base

After that:

.. code:: bash
    \pypsa-earth   $ conda install -c conda-forge mamba

Finally:

.. code:: bash
    \pypsa-earth   $ mamba env create -f envs/environment.yaml

Now everything should be ready to run the tutorial.

Run the tutorial model
---------------------

A tutorial data kit was developed to facilitate exploring the model.
The user can explore the majority of the model's functions on a local machine by running the tutorial, which uses fewer computational resources than the entire model does. 
Currently, the tutorial case study refers to a microgrid in Nigeria whose coverage area is defined by a rectangle whose vertices have the following cooridnates:

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
The tutorial model  run simulation will take a while (about 20..50 minutes).


Run a different case study
--------------------------

In this section you will find a small guide to learn how to configure the tool for your specific case study. 
As an example, we will refer to the study of a microgrid in Kenya whose coverage area is defined by a rectangle with vertices in these coordinates: 

- lon_max: 41.1141
- lon_min: 41.1086
- lat_min: -2.0596
- lat_max: -2.0526
To better understand this example, it might be helpful to read the section on configuration in the config.yaml file.

.. note::
    To find the coordinates of a specific study area, one functional way is to take advantage of the web application OpenStreetMap.

First of all, it is always a good idea to start by running the default tutorial case study; this will allow you to see if the tool is working properly. 
Furthermore, following the run of the tutorial case, you will see a "config.yaml" file appear in the folder. 
This file conjointly with the config.pypsa-earth.yaml allows the user to specify the scenario they wish to analyze before running the algorithm.

The following commands are the most relevant for setting up your personal case study:

1. Configure the country: in the two configuration files ( config.yaml and config.pypsa-earth.yaml) you must enter under "country" the abbreviation of the country where your case study is located.In our exemple "KE"

.. code:: yaml

    countries: ["KE"]

.. note::
    If the nation you have chosen is not among the nations that are enabled for tutorial configuration ( currently the only ones available are Nigeria (NG), Benin (BJ) , Botswana (BW) and Morocco (MA) ) you must:
    - Under tutorial ( at the top of the code ): replace true with false
    - Replace in both configuration files "cutout-2013-era5-tutorial" with "cutout-2013-era5" ( you can easily do this by exploiting the replace command in the file [ctrl+H])


2. Configure enable section : this section of the file ensure the download of essential open-source data, including databundle and cost data. 
   In our case it is convenient to go to set: 
    - build_cutout = false
    - build_nature_raster = false
   
   In particular the built_cutout = false is due to the fact that ,for the wheater year ,2013 will be left and in this case the precompiled cutouts are automatically downloaded with the retrive_databoundle rule. 
   The choice of build_nature_raster = false, is due to the fact that in this way a precompiled file "data/nature.tiff",also downloaded with the databundle, is used.

   .. code:: yaml

    enable:
        retrieve_databundle: true  #  Recommended 'true', for the first run. Otherwise data might be missing.
        retrieve_cost_data: true  
        download_osm_data: true 
        build_cutout: false
        build_natura_raster: false 
        

3. Enter the coordinates of the microgrid: in the config.yaml file, in the  microgrid_list section you have to insert the microgrid coordinates. In our case, the coordinates are:

   .. code:: yaml

    microgrids_list:
    microgrid_1:
        lon_max: 41.1141
        lon_min: 41.1086
        lat_min: -2.0596
        lat_max: -2.0526

These are the basic commands for configuring another case study than the one selected in the default configuration file. 
Of course, there are many other items to further customize your case study. 
If you would like to explore this further, we recommend you take a look at the Pypsa-Earth documentation in the configuration section.

At this point we are ready to run the newly configured case study.
As in the previous case, to run the code use this command:

.. code:: bash
    .../pypsa-distribution (pypsa-earth) % snakemake -j 1 solve_network

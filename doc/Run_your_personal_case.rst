.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _Run_your_personal_case:

#######################
Run_your_personal_case
#######################

After reading the previous sections, you are able to customise the configuration file and are almost ready to launch your modelling. 
Only a few small operations remain to be done. This section will explain the latest updates and how to run.
To better explain what needs to be done, we have chosen to go through the various operations by referring to an example.

Suppose we want to run the simulation for a small region of Turkey, and 
suppose we want to run a simulation on an annual time horizon using the RAMP demand simulation tool.

.. note::

    Remember: when using the tool for the first time, it is always a good idea to run the tutorial case preset in the configuration file. 
    This allows you to see if the code is working correctly and generates the config.yaml file which you can edit to run your own simulation.

Choose the type of simulation:
--------------------------------------

In our case, the goal is to perform a realistic simulation and for an area outside those that can be simulated with the tutorial mode. 
For this, the first thing to do is to set the tutorial mode false.
So, let's open the config.yaml file and start with the changes:

.. code:: yaml

    tutorial: false


Specify the country of interest:
--------------------------------------
To select the region of interest, simply change the ISO 3166-1 alpha-2 code, 
which represents the official country abbreviations. In the case of Turkey: TR.

.. code:: yaml

    countries: ["TR"]


Configure `coordinates`
--------------------------
One of the most important tasks is certainly to correctly select the location of the microgrid. 

.. note::
    A useful strategy for knowing the coordinates to be entered in the code is to use www.openstreetmap.org. 
    Once you have selected the area of interest by moving around the map, press the ‘export’ button. 
    You should now see a window on the left with a set of co-ordinates. By selecting ‘Select a different area manually’ 
    you can change the boundaries of the micro-network and find the desired set of co-ordinates to insert into the tool.

To enter the coordinates, simply change the values in the microgrid_list section to the desired values, 
in our case:

.. code:: yaml

    microgrids_list:
        your_microgrid:
            lon_max: 27.29261
            lon_min: 27.23579
            lat_min: 37.81941
            lat_max: 37.88786


The code is able to simulate several systems at the same time. 
To do so, simply add an additional analogous block to your_microgrid. For example:
.. code:: yaml

    microgrids_list:
        microgrid_1:
            lon_max: 27.29261
            lon_min: 27.23579
            lat_min: 37.81941
            lat_max: 37.88786

        microgrid_2:
            lon_max: 27.3941
            lon_min: 27.3556
            lat_min: 37.9334
            lat_max: 37.9611

Configure ``enable`` section to download/build data
---------------------------------------------------------
This section determines how much of the data required for the simulation is to be obtained. 
A configuration of this type is effective for launching the simulation, the inputs can of course be changed as required. 
Reference can also be made to the documentation in PyPSA-Earth on this.

.. code:: yaml

    enable:
        retrieve_databundle: true  # Recommended 'true', for the first run. Otherwise data might be missing.
        retrieve_cost_data: true   # If true, it retrieves cost data from technology data and saves in resources/costs.csv, if false uses cost data in data/costs.csv
        download_osm_data: true  # If true, OpenStreetMap data will be downloaded for the selected countries
        download_osm_buildings: true  # If true, OpenStreetMap buildings will be downloaded for the selected countries
        download_osm_method: overpass # or earth_osm
        build_cutout: false
        build_natura_raster: false  # If True, than an exclusion raster will be build


Specify the weather year scope
------------------------------
Let us assume, as already mentioned, that we want to run a simulation on an annual horizon. 
With this operation, the snapshots that will be associated with the network are set.

.. note::
    Attention: with this operation, the year of the climate data used for the calculation of the renewable yield is not selected, 
    this information depends on the cutout used. 

.. code:: yaml

    snapshots:
        start: "2013-01-01"
        end: "2014-01-01"
        inclusive: "left" # end is not inclusive

specify the method of load calculation :
------------------------------
Currently, energy demand can be simulated with two different strategies. 
Let us assume in this case that we want to use the strategy based on RAMP. 
With this modelling, the utility is divided into five representative classes, the amount of population associated with each "tier" is determined with the config.yaml parameter "tier_percent".
Inside the data/ramp folder are the excel files with the representative parameters of each class, which can be modified if necessary.

.. code:: yaml

    build_demand_type:
        type: "Ramp"
        std: "on"


Updates of configurations required for simulations not in tutorial mode:
------------------------------
.. note::

   Pypsa-Distribution is a tool under development, we are constantly working to improve the usability of the tool! 
   Currently, to run simulations in different areas than in the tutorial, it is also necessary to make these small changes.
   This is an interim solution which we aim to make more functional as soon as possible. 

Currently, in order to be able to run the simulation correctly, it is also necessary to slightly modify the config.pypsa-earth.yaml file; 
we are working on making this part of the procedure no longer necessary, but for the time being it is:
The changes you need to make to the config.pypsa-earth.yaml file are:
- Change the type of simulation from tutorial true to tutorial false:
  .. code:: yaml

    tutorial: false
- Change the country of interest:
  .. code:: yaml

    countries: ["TR"]
- Replace each time cutout-2013-era-tutorial appears in the config.yaml file with cutout-2013-era. 
  If you prefer to use another cutout, replace the name of the cutout above with the one you are interested in.
  
Run the simulation:
------------------------------
You are now officially ready to run the simulation, remember to open a wsl terminal and use the command to execute the run!

.. code:: bash
    .../pypsa-distribution (pypsa-earth) % snakemake -j 1 solve_network

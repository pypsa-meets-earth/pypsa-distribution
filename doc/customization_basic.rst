.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _customization_basic:

#######################
Basic customization
#######################

The goal of this section is to show you how to make basic settings to customize the configuration file and consequently the execution of the tool.
First, it is worth spending a few words on the configuration file.
This type of file is used to define and customize the settings and parameters for the execution of the workflow or process, so any customization of the tool will be performed on that file. 
Within the work folder there is always more than one configuration file, Before the first run of any rule within the code, you will see two configuration files in the folder: 
- config.distirbution.yaml 
- config.pypsa-earth.yaml.
When you run the tool, a new file will be created in the folder: config.yaml.
This is the file that contains all the specifications that will be executed in your simulation,
and it is therefore the file that you must take care to edit correctly.

If this is the first time you have approached these matters, you may be wondering why such a complex structure, the reason being that the config.yaml file that is generated is not a file that is tracked 
differently from the other two. This makes it possible to work in such a way that the config.yaml file can be modified by the user whenever he wants to run a simulation, 
while the config.distirbution.yaml and config.pypsa-earth.yaml files are only modified when the user wants to propose structural changes to the code.


In the next steps we will see the basic settings in order to customize the configuration file, in case you are interested in more in-depth customization, you can refer to the configuration section in the Pypsa-Earth documentation.

Choose the type of simulation:
--------------------------------------
It is possible to run simulations with the Pypsa Distribution tool in two modes:

- tutorial mode: this is a faster simulation in which less data is used. It was created to allow users to become familiar with the tool, without necessarily having to run heavy simulations.
  Remember: For realistic case studies it will be necessary to run the tool in full configuration and not in tutorial one.
  
- full ( non-tutorial) mode: in this case a full simulation is performed and the use of committed resources is greater.

Since the tutorial mode is designed only to become familiar with the tool, and not to perform analysis in a broad sense, it is possible to run the simulation in tutorial mode located in these nations: 

- Nigeria [NG]
- Benin [BJ]
- Botswana [BW]
- Morocco [MA]

The tutorial mode is selectable by placing within the config file:

.. code:: yaml

    tutorial: true

To perform a full simulation it will be necessary to set:

.. code:: yaml

    tutorial: false


Specify the country of interest:
--------------------------------------

Selection of the country of interest it's possible using the "coutries" argument at the top of the config file.

.. code:: yaml

    countries: ["NG"]


Please remember,before proceeding, that as has already been mentioned, in case you are running a simulation in tutorial mode you will only be able to select one nation from NG,BJ,BW or MA; 
and that for realistic case studies it will be necessary to run the tool in full configuration

Configure `coordinates`
--------------------------
Currently, the tool allows us to select the coverage area of interest by going to enter the coordinates of the vertices of the rectangular surface we want to be covered.
At the end of the config.yaml file you will find the microgrid list section where you can go to enter the coordinates. 
You can take advantage of the webb OpenStreetMap application to generate the coordinates.
The code offers the possibility of analysing more than one area of interest at a time.

.. code:: yaml

    microgrids_list:
        your_microgrid:
            lon_max: 41.1141
            lon_min: 41.1086
            lat_min: -2.0596
            lat_max: -2.0526


Configure ``enable`` section to download/build data
---------------------------------------------------------
The tool bases its operation on specific data obtained from open-source platforms. 
For a successful model run, ensure the download of essential open-source data, including databundle and cost data, is activated in the ``enable`` section:

.. code:: yaml

    enable:
        retrieve_databundle: true  # Recommended 'true', for the first run. Otherwise data might be missing.
        retrieve_cost_data: true   # If true, it retrieves cost data from technology data and saves in resources/costs.csv, if false uses cost data in data/costs.csv
        download_osm_data: true  # If true, OpenStreetMap data will be downloaded for the selected countries
        download_osm_buildings: true  # If true, OpenStreetMap buildings will be downloaded for the selected countries
        download_osm_method: overpass # or earth_osm
        build_cutout: false
        build_natura_raster: false  # If True, than an exclusion raster will be build
        
        # If "build_cutout" : true, then environmental data is extracted according to `snapshots` date range and `countries`
          # prepare_links_p_nom: false


Regarding topics in the enable section:

- Retrive_databundle and retrive_cost_data: After the initial run, it is a good idea to set retrive databundle and cost data to false to avoid re-downloading unnecessary data.
- Dowload_osm_data and download_osm_buildings: These two options are used to download OpenStreetMap data for the selected countries. The download_osm_method option allows you to choose the method to use to download the data. The overpass method is the default method, and currently the most effective.
- Build_cutout: When using weather year 2013 it is advisable to set "build_cutout: false" because the precompiled cutouts are automatically downloaded with the "retrive_databundle: true". When using a weather year other than 2013 it is essential to set "build_cutout: true" to generate custom cutouts.
    Caution: when using the Build_cutout rule, it is essential to first configure the Copernicus Climate Data Store API ( read the instructions).
    After the first run and successful custom cutout generation, build_cutout can be switched to false to avoid rebuilding the cutout.
- Build_nature_raster: When "build_nature_raster" is configured to false, the exclusion raster for protected areas is taken from the precompiled file "data/nature.tiff" downloaded with the databundle.Conversely, if "build_nature_raster" is set to "true" the exclusion raster is calculated using the "build_nature_raster" rule.
  After the initial run, it is recommended to set the retrieval of databundle and cost data to ``false`` to prevent unnecessary redownloading of data.


specify simulation time horizon:
----------------------------------
With these arguments, it is possible to set the time horizon of the simulation
Likewise, the example's temporal scope can be restricted (e.g. to 7 days):

.. code:: yaml

    snapshots:
        start: "2013-03-01"
        end: "2013-03-07"
        inclusive: "left" # end is not inclusive

.. note::

    Ensure that the selected date range aligns with the dates available in the cutout dataset. If the weather data within the cutouts corresponds to the year 2013, 
    then the range of snapshots should fall within that same year.
    In in the configuration file you will also find the variable "year" this specifies the reference year for the local population data. 
    It is your choice whether or not to keep the data temporally concordant

specify the method of load calculation :
-----------------------

This section specifies some parameters needed to generate demand profiles. 

.. code:: yaml

    build_demand_type:
        type: "From_file"
        std: "on"

This makes it possible to select the methodology adopted by the code for estimating demand. 
When "From_file" is selected, a preset consumption profile is used, which is scaled appropriately for the inhabitants of the area. 
Selecting "Ramp" makes it possible to use the RAMP tool for the generation of costumised usage profiles. 
For more in-depth information, it may be useful to see the information page: https://rampdemand.org/

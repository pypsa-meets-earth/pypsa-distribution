.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _customization_basic1:

#######################
2. Basic customization
#######################

The goal of this section is to show you how to make basic settings to customize the configuration file and consequently the execution of the tool.
First, it is worth spending a few words on the configuration file.
Within the work folder there is always a configuration file ( typically the file that is read by snakemake is the "config.yaml).
This type of file is used to define and customize the settings and parameters for the execution of the workflow or process, so any customization of the tool will be performed on that file. 
In the next steps we will see the basic settings in order to customize the configuration file, in case you are interested in more in-depth customization, you can refer to the configuration section in the Pypsa-Earth documentation.


Specify the country/region of interest
--------------------------------------

The model can be adapted to include any country, multiple countries (e.g. `Nigeria` and `Benin`)  using ``countries`` argument:

.. code:: yaml

    countries: ["NG"]

Configure ``enable`` section to download/build data
---------------------------------------------------------
The tool bases its operation on specific data obtained from open-source platforms. 
For a successful model run, ensure the download of essential open-source data, including databundle and cost data, is activated in the ``enable`` section:

.. code:: yaml

    enable:
        retrieve_databundle: true  #  Recommended 'true', for the first run. Otherwise data might be missing.
        retrieve_cost_data: true  # true: retrieves cost data from technology data and saves in resources/costs.csv, false: uses cost data in data/costs.csv
        download_osm_data: true  # If 'true', OpenStreetMap data will be downloaded for the above given countries
        build_cutout: false
        build_natura_raster: false # If True, than an exclusion raster will be build
        
        # If "build_cutout" : true, then environmental data is extracted according to `snapshots` date range and `countries`

Regarding topics in the enable section:

- Retrive_databundle and retrive_cost_data: After the initial run, it is a good idea to set retrive databundle and cost data to false to avoid re-downloading unnecessary data.
- Build_cutout: When using weather year 2013 it is advisable to set "build_cutout: false" because the precompiled cutouts are automatically downloaded with the "retrive_databundle: true". When using a weather year other than 2013 it is essential to set "build_cutout: true" to generate custom cutouts.
    Caution: when using the Build_cutout rule, it is essential to first configure the Copernicus Climate Data Store API ( read the instructions).
    After the first run and successful custom cutout generation, build_cutout can be switched to false to avoid rebuilding the cutout.
- Build_nature_raster: When "build_nature_raster" is configured to false, the exclusion raster for protected areas is taken from the precompiled file "data/nature.tiff" downloaded with the databundle.Conversely, if "build_nature_raster" is set to "true" the exclusion raster is calculated using the "build_nature_raster" rule.
  After the initial run, it is recommended to set the retrieval of databundle and cost data to ``false`` to prevent unnecessary redownloading of data.


Specify the weather year scope
------------------------------
With these arguments, it is possible to set the time horizon of the simulation
Likewise, the example's temporal scope can be restricted (e.g. to 7 days):

.. code:: yaml

    snapshots:
        start: "2013-03-01"
        end: "2013-03-07"
        inclusive: "left" # end is not inclusive

.. note::

    Ensure that the selected date range aligns with the dates available in the cutout dataset. If the weather data within the cutouts corresponds to the year 2013, then the range of snapshots should fall within that same year.

Specify the demand year
-----------------------

This section specifies some parameters needed to generate demand profiles. 

.. code:: yaml

    load_options:
      ssp: "ssp2-2.6"
      weather_year: 2013
      prediction_year: 2030
      scale: 1

The arguments you see are relative:
- Weather_year: sets the year referenced by the weather data for calculating electricity demand profiles for the selected area.
- Prediction_year: sets the Shared Socioeconomic Pathways (SSP) trajectory. Pypsa Earth uses the SSP2-2.6 scenario characterized by average challenges for mitigation and adaptation efforts to avoid global resgliding of about 2.6° by the end of the 21st century. Available values for weather_year and prediction_year can be checked by consulting the pypsa-earth/data/ssp2-2.6 folder
	Currently, pre-calculated demand data are available for weather years 2011, 2013, 2018 and prediction years 2030, 2040, 2050 and 2100.

Configure `atlite` section
--------------------------

To accurately model both temporally and spatially renewable availabilities such as wind and solar energy, historical climate data are processed with the atlite package.

.. code:: yaml

    atlite:
        nprocesses: 4
        cutouts:
            cutout-2013-era5:
                module: era5
                dx: 0.3  # cutout resolution
                dy: 0.3  # cutout resolution
                # The cutout time is automatically set by the snapshot range.

When you use precompiled cutouts, no editing of this section is required. 
However, when using precompiled cutouts, you must replace all "cutout-2013-era5" entries with the name of the custom cutout.
E.g.: if you simulate Kazakhstan with cutout: asia-2013-era5, each occurrence of cutout-2013-era5 should be updated to asia-2013-era5, which refers to the asia-2013-era5.nc file generated in the cutout folder.

Configure `coordinates`
--------------------------
Currently, the tool allows us to select the microgrid coverage area by going to enter the coordinates of the vertices of the rectangular surface we want to be covered by the microgrid.
At the end of the config.yaml file you will find the microgrid list section where you can go to enter the coordinates. 
You can take advantage of the webb OpenStreetMap application to generate the coordinates.

.. code:: yaml

    microgrids_list:
     microgrid_1:
        lon_max: 41.1141
        lon_min: 41.1086
        lat_min: -2.0596
        lat_max: -2.0526
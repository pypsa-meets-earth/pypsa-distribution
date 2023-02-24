# -*- coding: utf-8 -*-
"""
Estimates the population and the electric load of the microgrid.

Relevant Settings
-----------------
.. code:: yaml

    microgrids_list:
        Location:
            Centre:
                lon:
                lat:
            Sides:
                Deltalon:
                Deltalat:
        micA: 
            name:

    load:
        scaling_factor:

Inputs
------
- ``data/sample_profile.csv``: a load profile, which will be scaled through a scaling_factor to obtain the per person load

Outputs
-------
- ``data/Worldpop/population_file.tif: a tif file of the population of the selected country,
- ``resources/shapes/microgrid_shapes.geojson: a geojson file of the shape of the microgrid,
- ``resources/file_dir/country_masked.tif,
- ``resources/demand/microgrid_load.csv: the electric load of the microgid,
-------

Description
-----------

The rule :mod:`build_demand` contains functions that are used to create a shape file of the microgrid, to mask a raster with the shape file and to estimate 
the population. Then the population is multiplied for the per person load and the microgrid load is then obtained.

"""

import glob
import logging
import os
import shutil

import fiona
import geopandas as gpd
import pandas as pd
import rasterio
import rasterio.mask
import requests
from shapely.geometry import Polygon
from _helpers import configure_logging, get_country, sets_path_to_root
import json
_logger = logging.getLogger(__name__)
_logger.setLevel(logging.INFO)


def create_microgrid_shapes(microgrids_list, output_path): 

    """
    This function creates a rectangular shape of the microgrid and saves it as a .geojson file.
    The shape is defined by the coordinates of the angles of the rectangle.
    The resulting file is saved to the specified output_path.
    """
    
    microgrids_list=microgrids_list
    microgrids_list_df=pd.DataFrame(microgrids_list)

    microgrid_shapes = []
    microgrid_names = []

    for col in range(len(microgrids_list_df.columns)):

        values=microgrids_list_df.iloc[:, col]

        x1 = values[1] #lon_max
        y1 = values[3] #lat_max

        x2 = values[0] #lon_min
        y2 = values[3] #lat_max

        x3 = values[0] #lon_min
        y3 = values[2] #lat_min 

        x4 = values[1] #lon_max
        y4 = values[2] #lat_min


        microgrid_shape = Polygon([(x1, y1), (x2, y2), (x3, y3), (x4, y4), (x1, y1)])

        microgrid_name=f"microgrid_{col+1}"
        microgrid_shapes.append(microgrid_shape)
        microgrid_names.append(microgrid_name)

    microgrid_gdf = gpd.GeoDataFrame({'name': microgrid_names, 'geometry': microgrid_shapes})

    output_dict = json.loads(microgrid_gdf.to_json())
    output_json = json.dumps(output_dict, indent=4)

    with open(output_path, "w") as f:
        f.write(output_json)


def create_masked_file(raster_path, geojson_path, output_path):
    """
    This function masks a raster with a shape defined in a geojson file.
    The raster file is specified with the raster_path, the shape is specified with the geojson_path,
    and the resulting masked raster is saved to the specified output_path.
    """

    # Read the geojson file and convert it to a shapefile
    gdf = gpd.read_file(geojson_path)
    gdf.to_file("microgrid_shapes.shp")

    # Open the shapefile and extract the shape geometry
    with fiona.open("microgrid_shapes.shp", "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]

    # Open the raster and mask it using the shapes
    with rasterio.open(raster_path) as src:
        out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
        out_meta = src.meta

    # update the metadata for the output raster
    out_meta.update(
        {
            "driver": "GTiff",
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform,
        }
    )

    # Save the masked raster to the specified output path
    with rasterio.open(output_path, "w", **out_meta) as dest:
        dest.write(out_image)


def estimate_microgrid_population(masked_file, p, sample_profile, output_file):
    """
    This function estimates the population of a microgrid based on a mask file and a sample profile of electricity demand.
    The mask file is specified with the masked_file, the sample profile is specified with the sample_profile,
    and the output file for the estimated microgrid population is specified with the output_file.
    """

    with rasterio.open(masked_file) as fp:
        data = fp.read(1)
        pop_microgrid = data[data >= 0].sum()
        
    # Read the sample profile of electricity demand
    total_load = pd.read_csv(sample_profile)
    total_load = total_load["0"]

    # Calculate the per-person electricity demand and convert it as a pandas dataframe
    per_person_load = total_load * (1 / p)
    per_person_load = pd.DataFrame(per_person_load)

    # Calculate the microgrid electric load
    microgrid_load = per_person_load * pop_microgrid

    # Save the microgrid load to the specified output file
    microgrid_load.to_csv(output_file, index=False)

    return microgrid_load


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("build_demand")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)
    
    sample_profile = snakemake.input["sample_profile"]

    create_microgrid_shapes(
        snakemake.config["microgrids_list"], 
        snakemake.output["microgrid_shapes"],)

    create_masked_file(
        "data/Worldpop/population_file.tif",
        snakemake.output["microgrid_shapes"],
        snakemake.output["country_masked"],
    )

    estimate_microgrid_population(
        snakemake.output["country_masked"],
        snakemake.config["load"]["scaling_factor"],
        sample_profile,
        snakemake.output["electric_load"],
    )

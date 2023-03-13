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
- ``resources/shapes/microgrid_shape.geojson: a geojson file of the shape of the microgrid,
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
from _helpers_dist import (
    configure_logging,
    sets_path_to_root,
    two_2_three_digits_country,
)

_logger = logging.getLogger(__name__)
_logger.setLevel(logging.INFO)


def create_microgrid_shape(xcenter, ycenter, DeltaX, DeltaY, name, output_path):
    """
    This function creates a rectangular shape and saves it as a .geojson file.
    The shape is defined by its center coordinates (xcenter, ycenter) and its dimensions (DeltaX, DeltaY).
    The shape is also given a name, which is saved as a property in the .geojson file.
    The resulting file is saved to the specified output_path.
    """

    # Define the coordinates of the corners of the rectangle
    x1 = xcenter - DeltaX * 0.5
    y1 = ycenter + DeltaY * 0.5

    x2 = xcenter + DeltaX * 0.5
    y2 = ycenter + DeltaY * 0.5

    x3 = xcenter - DeltaX * 0.5
    y3 = ycenter - DeltaY * 0.5

    x4 = xcenter + DeltaX * 0.5
    y4 = ycenter - DeltaY * 0.5

    microgrid_name = name

    my_feature = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {"name": "microgrid_name"},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[[x1, y1], [x2, y2], [x3, y3], [x4, y4]]],
                },
            },
        ],
    }

    # Convert the feature to a GeoDataFrame
    gdf = gpd.GeoDataFrame.from_features(my_feature)

    # Save the GeoDataFrame to a .geojson file
    gdf.to_file(output_path)


def get_WorldPop_path(
    country_code,
    year,
    out_logging,
):
    """
    Download tiff file for each country code using the standard method from worldpop datastore with 1kmx1km resolution.

    Parameters
    ----------
    country_code : str
        Two letter country codes of the downloaded files.
        Files downloaded from https://data.worldpop.org/ datasets WorldPop UN adjusted
    year : int
        Year of the data to download
    Returns
    -------
    WorldPop_inputfile : str
        Path of the file
    """

    if out_logging:
        _logger.info("Download WorldPop datasets")

    three_digits_code = two_2_three_digits_country(country_code)

    return os.path.join(
        os.getcwd(),
        "pypsa-earth",
        "data",
        "Worldpop",
        f"{three_digits_code.lower()}_ppp_{year}_UNadj_constrained.tif",
    )  # Input filepath tif


def create_masked_file(raster_path, geojson_path, output_path):
    """
    This function masks a raster with a shape defined in a geojson file.
    The raster file is specified with the raster_path, the shape is specified with the geojson_path,
    and the resulting masked raster is saved to the specified output_path.
    """

    # Read the geojson file and convert it to a shapefile
    gdf = gpd.read_file(geojson_path)

    # Open the raster and mask it using the shapes
    with rasterio.open(raster_path) as src:
        out_image, out_transform = rasterio.mask.mask(src, gdf.geometry, crop=True)
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
        from _helpers_dist import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("build_demand")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)

    sample_profile = snakemake.input["sample_profile"]

    create_microgrid_shape(
        snakemake.config["microgrids_list"]["Location"]["Centre"]["lon"],
        snakemake.config["microgrids_list"]["Location"]["Centre"]["lat"],
        snakemake.config["microgrids_list"]["Location"]["Sides"]["Deltalon"],
        snakemake.config["microgrids_list"]["Location"]["Sides"]["Deltalat"],
        snakemake.config["microgrids_list"]["micA"]["name"],
        snakemake.output["microgrid_shape"],
    )

    worldpop_path = get_WorldPop_path(
        snakemake.config["countries"][0],  # TODO: this needs fix to generalize the countries
        snakemake.config["year"],
        False,
    )

    create_masked_file(
        worldpop_path,
        snakemake.output["microgrid_shape"],
        snakemake.output["country_masked"],
    )

    estimate_microgrid_population(
        snakemake.output["country_masked"],
        snakemake.config["load"]["scaling_factor"],
        sample_profile,
        snakemake.output["electric_load"],
    )

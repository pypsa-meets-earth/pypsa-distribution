# -*- coding: utf-8 -*-
"""
Estimates the population and the electric load of each microgrid.
Relevant Settings
-----------------
.. code:: yaml
    microgrids_list:
        microgridX: 
          lon_min:
          lon_max: 
          lat_min: 
          lat_max: 
    load:
        scaling_factor:
Inputs
------
- ``data/sample_profile.csv``: a load profile, which will be scaled through a scaling_factor to obtain the per person load
Outputs
-------
- ``resources/shapes/microgrid_shapes.geojson: a geojson file of the shape of each microgrid,
- ``resources/masked_files/masked_file_{i+1}.tif,
- ``resources/demand/microgrid_load_{i+1}.csv: the electric load of the microgid,
-------
Description
-----------
The rule :mod:`build_demand` contains functions that are used to create a shape file of the microgrid, to mask a raster with the shape file and to estimate 
the population. Then the population is multiplied for the per person load and the microgrid load is then obtained. The process applies to all the microgrids specified in config.yaml.
"""

import logging
import os

import geopandas as gpd
import pandas as pd
import pypsa
import rasterio
import rasterio.mask
from _helpers_dist import (
    configure_logging,
    sets_path_to_root,
    two_2_three_digits_country,
)

_logger = logging.getLogger(__name__)
_logger.setLevel(logging.INFO)


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
        "WorldPop",
        f"{three_digits_code.lower()}_ppp_{year}_UNadj_constrained.tif",
    )  # Input filepath tif


def estimate_microgrid_population(
    n, p, raster_path, shapes_path, sample_profile, output_file
):
    # Read the sample profile of electricity demand and extract the column corresponding to the electric load
    per_unit_load = pd.read_csv(sample_profile)["0"] / p

    # Dataframe of the load
    microgrid_load = pd.DataFrame()

    # Load the GeoJSON file with the shapes to mask the raster
    shapes = gpd.read_file(shapes_path)

    # Mask the raster with each shape and save each masked raster as a new file
    for i, shape in shapes.iterrows():
        with rasterio.open(raster_path) as src:
            # Mask the raster with the current shape
            masked, out_transform = rasterio.mask.mask(src, [shape.geometry], crop=True)
            out_meta = src.meta.copy()
            out_meta.update(
                {
                    "driver": "GTiff",
                    "height": masked.shape[1],
                    "width": masked.shape[2],
                    "transform": out_transform,
                }
            )

        pop_microgrid = masked[masked >= 0].sum()

        col_name = "new_bus_microgrid_" + str(i + 1)
        microgrid_load[col_name] = per_unit_load * pop_microgrid

    # Save the microgrid load to a CSV file with snapshots index
    microgrid_load.insert(0, "snapshots", n.snapshots)
    microgrid_load.set_index("snapshots", inplace=True)
    microgrid_load.to_csv(output_file, index=True)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_dist import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("build_demand")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.create_network)
    sample_profile = snakemake.input["sample_profile"]

    assert (
        len(snakemake.config["countries"]) == 1
    ), "Error: only a country shall be specified"

    worldpop_path = get_WorldPop_path(
        snakemake.config["countries"][
            0
        ],  # TODO: this needs fix to generalize the countries
        snakemake.config["year"],
        False,
    )

    estimate_microgrid_population(
        n,
        snakemake.config["load"]["scaling_factor"],
        worldpop_path,
        snakemake.input["microgrid_shapes"],
        sample_profile,
        snakemake.output["electric_load"],
    )

# -*- coding: utf-8 -*-

# TODO: Add docstring

import json
import os

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _helpers_dist import configure_logging, sets_path_to_root
from shapely.geometry import Point, Polygon


def extract_points(microgrid_shape_path, buildings_path, output_path):
    """
    From the downloaded data, extracts buildings located within the boundaries of each microgrid geometry 
    and associates them with the respective microgrid name.

    Parameters
    ----------
    microgrid_shape_path : str
        Path to the GeoJSON file containing microgrid geometries.
    buildings_path : str
        Path to the GeoJSON file containing building geometries.
    output_path : str
        Path where the resulting GeoJSON file will be saved.

    Returns
    -------
    GeoDataFrame
        A GeoDataFrame containing the filtered buildings with an added field "name_microgrid" 
        that associates each building to its corresponding microgrid.
    """
    
    # Load the GeoJSON files
    microgrid = gpd.read_file(microgrid_shape_path)
    buildings = gpd.read_file(buildings_path)
    # Create a GeoDataFrame to accumulate the results
    result = gpd.GeoDataFrame(columns=buildings.columns)
    # Iterate over each microgrid geometry
    for idx, microgrid_shape in microgrid.iterrows():
        # Extract the name of the microgrid
        microgrid_name = microgrid_shape["name"]
        # Filter buildings located within the microgrid geometry
        buildings_in_microgrid = buildings[buildings.geometry.within(microgrid_shape.geometry)]
        # Add or replace the "name_microgrid" field with the microgrid name
        buildings_in_microgrid = buildings_in_microgrid.copy()
        buildings_in_microgrid["name_microgrid"] = microgrid_name
        # Append the filtered buildings to the final result
        result = gpd.GeoDataFrame(pd.concat([result, buildings_in_microgrid], ignore_index=True))
    # Save the final result as a GeoJSON file
    result.to_file(output_path, driver="GeoJSON")

    return result


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_dist import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("clean_earth_osm_data")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)

    extract_points(
        snakemake.input["microgrid_shapes"],
        snakemake.input["all_buildings"],
        snakemake.output["microgrid_building"],
    )

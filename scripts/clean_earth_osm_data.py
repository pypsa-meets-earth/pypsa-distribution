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
    microgrid = gpd.read_file(microgrid_shape_path)
    xmin, ymin, xmax, ymax = microgrid.total_bounds

    buildings = gpd.read_file(buildings_path)
    buildings_in_microgrid = buildings.cx[xmin:xmax, ymin:ymax]

    buildings_in_microgrid.to_file(output_path)

    return buildings_in_microgrid


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

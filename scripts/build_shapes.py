# -*- coding: utf-8 -*-
import logging
import os

import geopandas as gpd
import pandas as pd
from _helpers_dist import configure_logging, save_to_geojson, sets_path_to_root
from shapely.geometry import Polygon

_logger = logging.getLogger(__name__)
_logger.setLevel(logging.INFO)


def create_microgrid_shapes(microgrids_list, output_path):
    """
    This function creates a rectangular shape of the microgrid and saves it as a .geojson file.
    The shape is defined by the coordinates of the angles of the rectangle.
    The resulting file is saved to the specified output_path.
    """

    microgrids_list = microgrids_list
    microgrids_list_df = pd.DataFrame(microgrids_list)

    microgrid_shapes = []
    microgrid_names = []

    for col in range(len(microgrids_list_df.columns)):
        values = microgrids_list_df.iloc[:, col]

        # Definition of the vertixes of the rectangle
        Top_left = (values[0], values[3])
        Top_right = (values[1], values[3])
        Bottom_right = (values[1], values[2])
        Bottom_left = (values[0], values[2])

        microgrid_shape = Polygon(
            [Top_left, Top_right, Bottom_right, Bottom_left, Top_left]
        )

        microgrid_name = f"microgrid_{col+1}"
        microgrid_shapes.append(microgrid_shape)
        microgrid_names.append(microgrid_name)

    microgrid_gdf = gpd.GeoDataFrame(
        {"name": microgrid_names, "geometry": microgrid_shapes}
    )

    save_to_geojson(microgrid_gdf, output_path)


def create_bus_regions(microgrids_list, output_path):
    """
    This function creates a geojson shape of the microgrid.
    The shape is defined by the coordinates of the angles of the rectangle and by another point,
    individuated by x and y which are the coordinates of the center of the microgrid.
    The resulting file is saved to the specified output_path.
    """

    microgrids_list = microgrids_list
    microgrids_list_df = pd.DataFrame(microgrids_list)

    microgrid_shapes = []
    microgrid_names = []
    microgrid_x = []
    microgrid_y = []

    for col in range(len(microgrids_list_df.columns)):
        values = microgrids_list_df.iloc[:, col]

        # Definition of the vertices of the rectangle
        Top_left = (values[0], values[3])
        Top_right = (values[1], values[3])
        Bottom_right = (values[1], values[2])
        Bottom_left = (values[0], values[2])

        microgrid_shape = Polygon(
            [Top_left, Top_right, Bottom_right, Bottom_left, Top_left]
        )

        # The bus is the central bus of each microgrid
        microgrid_name = f"new_bus_microgrid_{col+1}"
        microgrid_shapes.append(microgrid_shape)
        microgrid_names.append(microgrid_name)

        # Centre of the rectangle of the microgrid
        x = (values[0] + values[1]) * 0.5
        y = (values[2] + values[3]) * 0.5
        microgrid_x.append(x)
        microgrid_y.append(y)

    microgrid_gdf = gpd.GeoDataFrame(
        {
            "name": microgrid_names,
            "x": microgrid_x,
            "y": microgrid_y,
            "geometry": microgrid_shapes,
        }
    )

    save_to_geojson(microgrid_gdf, output_path)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_dist import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("build_shapes")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)

    create_microgrid_shapes(
        snakemake.config["microgrids_list"],
        snakemake.output["microgrid_shapes"],
    )

    create_bus_regions(
        snakemake.config["microgrids_list"],
        snakemake.output["microgrid_bus_shapes"],
    )

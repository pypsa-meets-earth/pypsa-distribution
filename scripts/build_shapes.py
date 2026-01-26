# -*- coding: utf-8 -*-
import logging
import os

import geopandas as gpd
import pandas as pd
from _helpers_dist import configure_logging, save_to_geojson, sets_path_to_root
from shapely.geometry import Polygon

_logger = logging.getLogger(__name__)
_logger.setLevel(logging.INFO)


def create_microgrid_shapes(microgrids_list, output_path, country_code):
    """
    Creates rectangular shapes for each microgrid in the list of microgrids in the config.yaml file
    and saves them as a GeoJSON file.
    Parameters
    ----------
    microgrids_list : dict
        Dictionary containing the microgrid names and their bounding box coordinates (lat_min, lon_min, lat_max, lon_max).

    output_path : str
       Path where the GeoJSON file will be saved.
    country_code : str
       Country code to store under the 'name' property.
    """

    # Open the input dictionary into a pandas DataFrame for easier processing
    microgrids_list_df = pd.DataFrame(microgrids_list)

    # Initialize lists to store shapes and names of each microgrid
    microgrid_shapes = []
    microgrid_names = []

    # Iterate over each column (representing a microgrid) in the DataFrame
    for col in range(len(microgrids_list_df.columns)):
        # Extract the bounds of the rectangle for the current microgrid
        values = microgrids_list_df.iloc[:, col]
        # Define the vertices of the rectangle
        Top_left = (values.iloc[0], values.iloc[3])
        Top_right = (values.iloc[1], values.iloc[3])
        Bottom_right = (values.iloc[1], values.iloc[2])
        Bottom_left = (values.iloc[0], values.iloc[2])
        # Create a Polygon shape from the rectangle's vertices
        microgrid_shape = Polygon(
            [Top_left, Top_right, Bottom_right, Bottom_left, Top_left]
        )
        # Assign a unique name to the microgrid based on its name in the config
        microgrid_name = f"microgrid_{col+1}"
        # Append the shape and name to the respective lists
        microgrid_shapes.append(microgrid_shape)
        microgrid_names.append(microgrid_name)

    # Create a GeoDataFrame with the collected names and shapes
    microgrid_gdf = gpd.GeoDataFrame(
        {"name_microgrid": microgrid_names, "geometry": microgrid_shapes}
    )
    microgrid_gdf["name"] = country_code

    # Save the GeoDataFrame to a GeoJSON file
    save_to_geojson(microgrid_gdf, output_path)


def create_bus_regions(microgrids_list, output_path, country_code):
    """
    Creates bus regions for each microgrid in the list of microgrids and saves them as a GeoJSON file.
    The generated shape will be used for the calculation of renewable energy producibility,
    which will be associated with the bus generated at the center of the geometry.
    Parameters
    ----------
    microgrids_list : dict
        Dictionary containing the microgrid names and their bounding box coordinates (lat_min, lon_min, lat_max, lon_max).

    output_path : str
       Path where the GeoJSON file will be saved.
    country_code : str
       Country code to store under the 'name' property.
    """

    # Open the input dictionary as pandas DataFrame for easier processing
    microgrids_list_df = pd.DataFrame(microgrids_list)

    # Initialize lists to store shapes, names, and coordinates
    microgrid_shapes = []
    microgrid_names = []
    microgrid_x = []  # Stores the x-coordinates of the centers of the rectangles
    microgrid_y = []  # Stores the y-coordinates of the centers of the rectangles

    # Iterate over each column in the DataFrame
    for col in range(len(microgrids_list_df.columns)):
        values = microgrids_list_df.iloc[:, col]
        microgrid_name = microgrids_list_df.columns[col] + "_gen_bus"

        # Define the vertices of the rectangle
        Top_left = (values.iloc[0], values.iloc[3])
        Top_right = (values.iloc[1], values.iloc[3])
        Bottom_right = (values.iloc[1], values.iloc[2])
        Bottom_left = (values.iloc[0], values.iloc[2])

        # Create a Polygon shape from the rectangle's vertices
        microgrid_shape = Polygon(
            [Top_left, Top_right, Bottom_right, Bottom_left, Top_left]
        )

        # Append the shape and name to the respective lists
        microgrid_shapes.append(microgrid_shape)
        microgrid_names.append(microgrid_name)

        # Calculate the center of the rectangle
        x = (values.iloc[0] + values.iloc[1]) / 2
        y = (values.iloc[2] + values.iloc[3]) / 2
        microgrid_x.append(x)  # Append the x-coordinate of the center
        microgrid_y.append(y)  # Append the y-coordinate of the center

    # Create a GeoDataFrame from the collected names, shapes, and coordinates
    microgrid_gdf = gpd.GeoDataFrame(
        {
            "name_microgrid": microgrid_names,  # microgrid-specific name
            "x": microgrid_x,  # x-coordinates of the centers
            "y": microgrid_y,  # y-coordinates of the centers
            "geometry": microgrid_shapes,  # Polygon shapes of the regions
        }
    )
    microgrid_gdf["name"] = country_code

    # Save the GeoDataFrame to a GeoJSON file
    save_to_geojson(microgrid_gdf, output_path)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_dist import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("build_shapes")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)
    country_code = snakemake.params["countries"]

    create_microgrid_shapes(
        snakemake.config["microgrids_list"],
        snakemake.output["microgrid_shapes"],
        country_code=country_code,
    )

    create_bus_regions(
        snakemake.config["microgrids_list"],
        snakemake.output["microgrid_bus_shapes"],
        country_code=country_code,
    )

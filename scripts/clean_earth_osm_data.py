# -*- coding: utf-8 -*-

# TODO: Add docstring

import json
import geojson
import os

import matplotlib.pyplot as plt
import numpy as np
from _helpers_dist import configure_logging, sets_path_to_root
from shapely.geometry import Point, Polygon


def combine_geojson_files(folder_input_path, output_file):
    
    features = []

    for filename in os.listdir(folder_input_path):
        if filename.endswith(".geojson"):
            file_path = os.path.join(folder_input_path, filename)
            if os.stat(file_path).st_size == 0:
                print(f"Skipping empty file: {filename}")
                continue
            with open(file_path) as f:
                data = json.load(f)
                for feature in data["features"]:
                    features.append(feature)

    feature_collection = {"type": "FeatureCollection", "features": features}

    with open(output_file, "w") as f:
        json.dump(feature_collection, f)


def extract_points_inside_microgrids(
    input_building_file, input_microgrid_file, output_file
):
    with open(input_building_file) as f:
        points_geojson = json.load(f)

    with open(input_microgrid_file) as f:
        rectangle_geojson = json.load(f)

    # Iterate over the rectangles and extract the points inside each one
    points_inside_rectangles = []
    for feature in rectangle_geojson["features"]:
        rectangle_coords = feature["geometry"]["coordinates"][0]
        rectangle_coords = [
            (lon, lat) for lat, lon in rectangle_coords
        ]  # Swap lat and lon
        rectangle_polygon = Polygon(rectangle_coords)

        # Extract the points inside the rectangle
        points_in_rectangle = []
        for point_feature in points_geojson["features"]:
            point_coords = point_feature["geometry"]["coordinates"]
            point_coords = (point_coords[1], point_coords[0])
            point = Point(point_coords)
            if rectangle_polygon.contains(point):
                # Add the microgrid_id property to the point feature
                point_feature["properties"]["microgrid_id"] = feature["properties"][
                    "name"
                ]
                points_in_rectangle.append(point_feature)

        # Add the points to the list of points inside all rectangles
        points_inside_rectangles.extend(points_in_rectangle)

    # Save the results to a new geojson file
    result_geojson = {"type": "FeatureCollection", "features": points_inside_rectangles}

    with open(output_file, "w") as f:
        json.dump(result_geojson, f)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_dist import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("clean_earth_osm_data")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)
      
    combine_geojson_files(snakemake.input["folder_input_path"], snakemake.output["buildings_geojson"])

    extract_points_inside_microgrids(
        snakemake.output["buildings_geojson"],
        snakemake.input["microgrid_shapes"],
        snakemake.output["microgrids_buildings"],
    )

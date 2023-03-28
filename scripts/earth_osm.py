# -*- coding: utf-8 -*-
import json
import os

import matplotlib.pyplot as plt
import numpy as np
from _helpers_dist import configure_logging, sets_path_to_root
from scipy.spatial import Delaunay
from shapely.geometry import Point, Polygon
import json
import geojson


def transform_json_to_geojson(input_file, output_file):
    """
    Transform a custom JSON file of building data to a GeoJSON file.

    Parameters:
        input_file (str): The path to the input JSON file.
        output_file (str): The path to the output GeoJSON file.

    Returns:
        None

    Example:
        >>> transform_json_to_geojson("building.json", "output_building.geojson")

    """
    # Load the custom JSON data
    with open(input_file) as f:
        data = json.load(f)

    # Create a list of features from the nodes
    features = []
    for node_id, node_data in data['Data']['Node'].items():
        lon, lat = node_data['lonlat']
        geometry = geojson.Point((lon, lat))
        feature = geojson.Feature(id=node_id, geometry=geometry, properties=node_data['tags'])
        features.append(feature)

    # Create a FeatureCollection from the features
    feature_collection = geojson.FeatureCollection(features)

    # Write the FeatureCollection to a GeoJSON file
    with open(output_file, 'w') as f:
        geojson.dump(feature_collection, f)


def extract_points_inside_microgrids(input_building_file, input_microgrid_file, output_file):

    with open(input_building_file) as f:
        points_geojson = json.load(f)

    with open(input_microgrid_file) as f:
        rectangle_geojson = json.load(f)

# Iterate over the rectangles and extract the points inside each one
    points_inside_rectangles = []
    for feature in rectangle_geojson['features']:
        rectangle_coords = feature['geometry']['coordinates'][0]
        rectangle_coords = [(lat, lon) for lon, lat in rectangle_coords] # Swap lat and lon
        rectangle_polygon = Polygon(rectangle_coords)
    
    # Extract the points inside the rectangle
        points_in_rectangle = []
        for point_feature in points_geojson['features']:
            point_coords = point_feature['geometry']['coordinates']
            point_coords = (point_coords[1], point_coords[0])
            point = Point(point_coords)
            if rectangle_polygon.contains(point):
               # Add the microgrid_id property to the point feature
                point_feature['properties']['microgrid_id'] = feature['properties']['name']
                points_in_rectangle.append(point_feature)
    
    # Add the points to the list of points inside all rectangles
        points_inside_rectangles.extend(points_in_rectangle)

# Save the results to a new geojson file
    result_geojson = {
        'type': 'FeatureCollection',
        'features': points_inside_rectangles
    }

    with open(output_file, 'w') as f:
        json.dump(result_geojson, f)


def build_plot_from_points(input_file, output_path):


    with open(input_file) as f:
        data = json.load(f)

    # extract point coordinates and microgrid names from GeoJSON data
    points = []
    microgrid_names = []
    for feature in data["features"]:
        coords = feature["geometry"]["coordinates"]
        microgrid_name = feature["properties"].get("microgrid_id")
        if microgrid_name is not None:
            points.append(coords)
            microgrid_names.append(microgrid_name)

    # create separate triangulations for each microgrid
    triangulations = {}
    for microgrid_name in set(microgrid_names):
        # select points of the current microgrid
        indices = [i for i, name in enumerate(microgrid_names) if name == microgrid_name]
        microgrid_points = np.array([points[i] for i in indices])

        # create Delaunay triangulation for the current microgrid
        triangulations[microgrid_name] = Delaunay(microgrid_points)


    # Plot the points and triangulations for each microgrid
    for microgrid_name, tri in triangulations.items():
        # select points of the current microgrid
        indices = [i for i, name in enumerate(microgrid_names) if name == microgrid_name]
        microgrid_points = np.array([points[i] for i in indices])

        # plot the points
        plt.plot(microgrid_points[:, 0], microgrid_points[:, 1], "ko")

        # plot the triangulation
        plt.triplot(microgrid_points[:, 0], microgrid_points[:, 1], tri.simplices)

        # add title
        plt.title(f"Triangulation for {microgrid_name}")

        # save the plot
        plt.savefig(f"{output_path}/{microgrid_name}.png")
        
        # close the figure
        plt.close()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_dist import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("earth_osm")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)

    transform_json_to_geojson(snakemake.input["buildings_json"], 
                              snakemake.output["buildings_geojson"])

    extract_points_inside_microgrids(snakemake.output["buildings_geojson"],
                                     snakemake.input["microgrid_shapes"],
                                     snakemake.output["microgrids_buildings"])

    build_plot_from_points(snakemake.output["microgrids_buildings"],
                           snakemake.output["plot_delaunay"])

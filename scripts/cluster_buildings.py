# -*- coding: utf-8 -*-
import json
import logging
import os
from collections import Counter

import geopandas as gpd
import numpy as np
import pandas as pd
from _helpers_dist import (
    configure_logging,
    sets_path_to_root,
    two_2_three_digits_country,
)
from shapely.geometry import Point, shape
from sklearn.cluster import KMeans

_logger = logging.getLogger(__name__)
_logger.setLevel(logging.INFO)


def extract_points(input_file, output_file):
    microgrid_buildings = gpd.read_file(input_file)
    microgrid_buildings.rename(columns={"tags.building": "tags_building"}, inplace=True)
    features = []
    for row in microgrid_buildings.itertuples():
        if row.Type == "area":
            features.append(
                {
                    "properties": {
                        "id": row.id,
                        "type": row.Type,
                        "tags_building": row.tags_building,
                    },
                    "geometry": row.geometry,
                }
            )
    buildings_geodataframe = gpd.GeoDataFrame.from_features(features)
    buildings_geodataframe.to_file(output_file)


def get_central_points_geojson(input_filepath, output_filepath, n_clusters):
    microgrid_buildings = gpd.read_file(input_filepath)
    centroids_building = [
        (row.geometry.centroid.x, row.geometry.centroid.y)
        for row in microgrid_buildings.itertuples()
    ]
    microgrid_buildings["centroid_coordinates"] = centroids_building

    centroids_building = np.array(centroids_building)
    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(centroids_building)
    centroids = kmeans.cluster_centers_

    central_points = []
    for i in range(kmeans.n_clusters):
        cluster_points = centroids_building[kmeans.labels_ == i]
        distances = np.linalg.norm(cluster_points - centroids[i], axis=1)
        central_point_idx = np.argmin(distances)
        central_points.append(cluster_points[central_point_idx])

    central_features = []
    for i, central_point in enumerate(central_points):
        central_features.append(
            {
                "geometry": Point(central_point),
                "cluster": i,
            }
        )
    central_features = gpd.GeoDataFrame(central_features)
    central_features.to_file(output_filepath)


def get_central_points_geojson_with_buildings(
    input_filepath, output_filepath, n_clusters
):
    cleaned_buildings = gpd.read_file(input_filepath)

    centroids_building = [
        (row.geometry.centroid.x, row.geometry.centroid.y)
        for row in cleaned_buildings.itertuples()
    ]
    centroids_building = np.array(centroids_building)
    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(centroids_building)
    centroids = kmeans.cluster_centers_

    central_points = []
    for i in range(kmeans.n_clusters):
        cluster_points = centroids_building[kmeans.labels_ == i]
        distances = np.linalg.norm(cluster_points - centroids[i], axis=1)
        central_point_idx = np.argmin(distances)
        central_points.append(cluster_points[central_point_idx])

    central_features = []
    for i, central_point in enumerate(central_points):
        central_features.append(
            {"geometry": Point(central_point), "cluster": i, "buildings": []}
        )
    central_gdf = gpd.GeoDataFrame(central_features)

    for i, label in enumerate(kmeans.labels_):
        feature = cleaned_buildings.iloc[i]
        if feature.geometry.type == "Polygon":
            cluster_id = label
            building_tag = feature["tags_building"]
            central_gdf.at[cluster_id, "buildings"].append(building_tag)

    central_gdf["buildings"] = central_gdf["buildings"].apply(json.dumps)
    central_gdf.to_file(output_filepath, driver="GeoJSON")


def get_number_type_buildings(input_filepath, output_filepath):
    cleaned_buildings = gpd.read_file(input_filepath)
    counts = []
    for row in cleaned_buildings.itertuples():
        building_tag = row.buildings
        building_tag = json.loads(building_tag.replace("'", '"'))
        building_tag = pd.Series(building_tag)
        count = building_tag.value_counts()
        counts.append(count)
    counts = pd.DataFrame(counts).fillna(0).astype(int)
    counts["cluster"] = cleaned_buildings["cluster"].values
    counts.set_index("cluster", inplace=True)
    counts.to_excel(output_filepath)
    print("fino a qui tutto bene")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_dist import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("cluster_buildings")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)

    extract_points(
        snakemake.input["buildings_geojson"],
        snakemake.output["cleaned_buildings_geojson"],
    )

    get_central_points_geojson(
        snakemake.output["cleaned_buildings_geojson"],
        snakemake.output["clusters"],
        snakemake.config["buildings"]["n_clusters"],
    )

    get_central_points_geojson_with_buildings(
        snakemake.output["cleaned_buildings_geojson"],
        snakemake.output["clusters_with_buildings"],
        snakemake.config["buildings"]["n_clusters"],
    )

    get_number_type_buildings(
        snakemake.output["clusters_with_buildings"],
        snakemake.output["number_buildings_type"],
    )

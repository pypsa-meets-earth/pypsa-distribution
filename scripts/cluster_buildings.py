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


def buildings_classification(input_file, crs, house_area_limit):
    microgrid_buildings = gpd.read_file(input_file)
    microgrid_buildings.rename(columns={"tags.building": "tags_building"}, inplace=True)
    features = []
    for row in microgrid_buildings.itertuples():
        if row.geometry.type == "Polygon":
            features.append(
                {
                    "properties": {
                        "id": row.id,
                        "tags_building": row.tags_building,
                    },
                    "geometry": row.geometry,
                }
            )
    buildings_geodataframe = gpd.GeoDataFrame.from_features(features)
    microgrid_buildings = microgrid_buildings.to_crs(crs)
    area_m2 = microgrid_buildings.geometry.area
    area_m2 = pd.Series(area_m2)
    buildings_geodataframe["area_m2"] = area_m2
    for index, row in buildings_geodataframe.iterrows():
        if row.tags_building == "yes":
            if row.area_m2 < house_area_limit:
                buildings_geodataframe.loc[index, "tags_building"] = "house"
            else:
                buildings_geodataframe.loc[index, "tags_building"] = "yes"
    return buildings_geodataframe
    # buildings_geodataframe.to_file(output_file)


def get_central_points_geojson(
    input_filepath, output_filepath, n_clusters, crs, house_area_limit
):
    microgrid_buildings = buildings_classification(
        input_filepath, crs, house_area_limit
    )
    # microgrid_buildings = gpd.read_file(input_filepath)
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
    input_filepath, output_filepath, n_clusters, crs, house_area_limit
):
    cleaned_buildings = buildings_classification(input_filepath, crs, house_area_limit)
    # cleaned_buildings = gpd.read_file(input_filepath)
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

    features = []
    for i, row in enumerate(cleaned_buildings.itertuples()):
        if row.geometry.type == "Polygon":
            cluster_id = kmeans.labels_[i]
            features.append(
                {
                    "properties": {
                        "id": row.Index,
                        "cluster_id": cluster_id,
                        "tags_building": row.tags_building,
                        "area_m2": row.area_m2,
                    },
                    "geometry": row.geometry,
                }
            )
    buildings_geodataframe = gpd.GeoDataFrame.from_features(features)
    buildings_geodataframe.to_file(output_filepath)


def get_number_type_buildings(input_filepath, output_filepath):
    buildings_geodataframe = gpd.read_file(input_filepath)

    grouped_buildings = buildings_geodataframe.groupby("cluster_id")
    clusters = np.sort(buildings_geodataframe["cluster_id"].unique())
    counts = []
    for cluster in clusters:
        cluster_buildings = pd.DataFrame(grouped_buildings.get_group(cluster))
        building_tag = cluster_buildings["tags_building"]
        building_tag = pd.Series(building_tag)
        count = building_tag.value_counts()
        counts.append(count)
    counts = pd.DataFrame(counts).fillna(0).astype(int)
    counts["cluster"] = clusters
    counts.set_index("cluster", inplace=True)
    counts.to_csv(output_filepath)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_dist import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("cluster_buildings")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)

    crs = snakemake.params.crs["area_crs"]
    house_area_limit = snakemake.params.house_area_limit["area_limit"]

    get_central_points_geojson(
        snakemake.input["buildings_geojson"],
        snakemake.output["clusters"],
        snakemake.config["buildings"]["n_clusters"],
        crs,
        house_area_limit,
    )

    get_central_points_geojson_with_buildings(
        snakemake.input["buildings_geojson"],
        snakemake.output["clusters_with_buildings"],
        snakemake.config["buildings"]["n_clusters"],
        crs,
        house_area_limit,
    )

    get_number_type_buildings(
        snakemake.output["clusters_with_buildings"],
        snakemake.output["number_buildings_type"],
    )

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
    """
    Filters the data contained in all_raw_building, selecting only Polygon elements,
    after which the plan area is calculated for each building with the specified coordinate system
    and adds the information to the geodataframe.
    """
    microgrid_buildings = gpd.read_file(input_file)
    microgrid_buildings.rename(columns={"tags.building": "tags_building"}, inplace=True)
    microgrid_buildings = microgrid_buildings.loc[
        microgrid_buildings.geometry.type != "Point"
    ]
    microgrid_buildings = microgrid_buildings.to_crs(crs)
    microgrid_buildings["area_m2"] = microgrid_buildings.geometry.area
    idxs_house = microgrid_buildings.query(
        "(tags_building == 'yes') and (area_m2 < @house_area_limit)"
    ).index
    microgrid_buildings.loc[idxs_house, "tags_building"] = "house"
    return microgrid_buildings


def get_central_points_geojson_with_buildings(
    input_filepath,
    output_filepath_centroids,
    n_clusters,
    crs,
    house_area_limit,
    output_filepath_buildings,
    output_path_csv,
):
    """
    Divides the buildings into the desired number of clusters by using the kmeans function
    and returns three different output: a geodataframe with the coordinates of the centroids of each cluster,
    a dataframe with all of the buildings divided into clusters,
    a csv file where for each cluster the building types are counted
    """
    microgrid_buildings = buildings_classification(
        input_filepath, crs, house_area_limit
    )
    centroids_building = [
        (row.geometry.centroid.x, row.geometry.centroid.y)
        for row in microgrid_buildings.itertuples()
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
            {
                "geometry": Point(central_point),
                "cluster": i,
            }
        )
    central_features = gpd.GeoDataFrame(central_features, crs="EPSG:4326")
    central_features.to_file(output_filepath_centroids, driver="GeoJSON")

    clusters = []
    for i, row in enumerate(microgrid_buildings.itertuples()):
        cluster_id = kmeans.labels_[i]
        clusters.append(cluster_id)

    microgrid_buildings["cluster_id"] = clusters

    microgrid_buildings_gdf = gpd.GeoDataFrame(
        microgrid_buildings, crs=microgrid_buildings.crs
    )
    microgrid_buildings_gdf.to_file(output_filepath_buildings)
    buildings_class = pd.DataFrame(
        microgrid_buildings_gdf.groupby("cluster_id").tags_building.value_counts()
    )
    buildings_class.to_csv(output_path_csv)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_dist import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("cluster_buildings")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)

    crs = snakemake.params.crs["area_crs"]
    house_area_limit = snakemake.params.house_area_limit["area_limit"]

    get_central_points_geojson_with_buildings(
        snakemake.input["buildings_geojson"],
        snakemake.output["clusters"],
        snakemake.config["buildings"]["n_clusters"],
        crs,
        house_area_limit,
        snakemake.output["clusters_with_buildings"],
        snakemake.output["buildings_type"],
    )

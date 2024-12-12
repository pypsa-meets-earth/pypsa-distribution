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


def buildings_classification(input_file, crs):
    """
    Filters the data contained in all_raw_building, selecting only Polygon elements,
    after which the plan area is calculated for each building with the specified coordinate system
    and adds the information to the geodataframe.
    """
    microgrid_buildings = gpd.read_file(input_file)
    microgrid_buildings.rename(columns={"building": "tags_building"}, inplace=True)
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
    microgrids_list
):
    """
    Divides the buildings into the desired number of clusters by using the kmeans function
    and generates three outputs:
    - GeoJSON with the coordinates of the centroids of each cluster,
    - GeoJSON with all the buildings divided into clusters,
    - CSV file where the building types are counted for each cluster.
    """
    
    microgrid_buildings = buildings_classification(input_filepath, crs)

    
    all_central_features = gpd.GeoDataFrame(columns=["geometry", "cluster", "name_microgrid"])
    all_microgrid_buildings = gpd.GeoDataFrame(columns=microgrid_buildings.columns)
    all_buildings_class = pd.DataFrame()

    
    for grid_name, grid_data in microgrids_list.items():
        
        filtered_buildings = microgrid_buildings[microgrid_buildings["name_microgrid"] == grid_name]

    
        centroids_building = [
            (row.geometry.centroid.x, row.geometry.centroid.y)
            for row in filtered_buildings.itertuples()
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
                    "name_microgrid": grid_name,
                }
            )
        central_features_gdf = gpd.GeoDataFrame(
            central_features, crs=filtered_buildings.crs
        ).to_crs("EPSG:4326")
        all_central_features = pd.concat([all_central_features, central_features_gdf], ignore_index=True)

        
        clusters = kmeans.labels_
        filtered_buildings["cluster_id"] = clusters.astype(int)
        all_microgrid_buildings = pd.concat([all_microgrid_buildings, filtered_buildings], ignore_index=True)

       
        buildings_class = (
            filtered_buildings.groupby("cluster_id").tags_building.value_counts().reset_index(name="count")
        )
        buildings_class["name_microgrid"] = grid_name
        all_buildings_class = pd.concat([all_buildings_class, buildings_class], ignore_index=True)

    
    all_central_features.to_file(output_filepath_centroids, driver="GeoJSON")
    all_microgrid_buildings.to_file(output_filepath_buildings, driver="GeoJSON")
    all_buildings_class.to_csv(output_path_csv, index=False)



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
        snakemake.config["microgrids_list"],
    )

# -*- coding: utf-8 -*-
import json
import logging
import os
from collections import Counter
from pathlib import Path

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
    Filters the data contained in the input GeoJSON file, selecting only Polygon elements.
    Calculates the plan area for each building based on the specified coordinate system (CRS)
    and adds this information as a new column to the GeoDataFrame.
    Buildings classified as "yes" with an area below a predefined limit are reclassified as "house".

    Parameters
    ----------
    input_file : str
        Path to the input GeoJSON file containing building data.
    crs : str
        The coordinate reference system (CRS) to be used for area calculation.
    Returns
    -------
    GeoDataFrame
        A GeoDataFrame containing filtered and classified building data with the added "area_m2" column.
    """
    # Load the GeoJSON file
    microgrid_buildings = gpd.read_file(input_file)
    microgrid_buildings.rename(columns={"building": "tags_building"}, inplace=True)
    # Filter out elements that are Points, keeping only Polygons and MultiPolygons
    microgrid_buildings = microgrid_buildings.loc[
        microgrid_buildings.geometry.type != "Point"
    ]
    # Convert the GeoDataFrame to the specified CRS
    microgrid_buildings = microgrid_buildings.to_crs(crs)
    # Calculate the area of each building and store it in a new column "area_m2"
    microgrid_buildings["area_m2"] = microgrid_buildings.geometry.area
    # Identify buildings with "tags_building" = "yes" and area below the house_area_limit and reclassify these buildings as "house"
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
    microgrids_list,
    geo_crs
):
    """
    Divides buildings into a specified number of clusters using the KMeans algorithm and generates:
    - GeoJSON file containing the centroids of each cluster,
    - GeoJSON file containing all buildings with cluster assignments,
    - CSV file summarizing the count of building types within each cluster.

    Parameters
    ----------
    input_filepath : str
        Path to the input GeoJSON file containing building data.
    output_filepath_centroids : str
        Path to the output GeoJSON file for cluster centroids.
    n_clusters : int
        Number of clusters to divide the buildings into.
    crs : str
        The coordinate reference system (CRS) for spatial operations.
    house_area_limit : float
        The maximum area (in square meters) to classify a building as a "house".
    output_filepath_buildings : str
        Path to the output GeoJSON file containing clustered buildings.
    output_path_csv : str
        Path to the output CSV file summarizing building types per cluster.
    microgrids_list : dict
        Dictionary of microgrids with their names and bounding coordinates.

    """
    # Classify and process the buildings
    microgrid_buildings = buildings_classification(input_filepath, crs)

    # Prepare accumulators
    all_central_features = gpd.GeoDataFrame(
    columns=["geometry", "cluster", "name_microgrid"], crs=geo_crs
    )
    all_microgrid_buildings = gpd.GeoDataFrame(columns=microgrid_buildings.columns)
    all_buildings_class = pd.DataFrame()

    # Process each microgrid
    for grid_name, grid_data in microgrids_list.items():
        filtered_buildings = microgrid_buildings[
            microgrid_buildings["name_microgrid"] == grid_name
        ]

        # Extract building centroids
        centroids_building = [
            (row.geometry.centroid.x, row.geometry.centroid.y)
            for row in filtered_buildings.itertuples()
        ]
        centroids_building = np.array(centroids_building)

        # KMeans
        if isinstance(n_clusters, dict):
            kmeans = KMeans(n_clusters=n_clusters[grid_name], random_state=0).fit(
                centroids_building
            )
        else:
            kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(
                centroids_building
            )

        # Central points of clusters
        centroids = kmeans.cluster_centers_
        central_points = []
        for i in range(kmeans.n_clusters):
            cluster_points = centroids_building[kmeans.labels_ == i]
            distances = np.linalg.norm(cluster_points - centroids[i], axis=1)
            central_point_idx = np.argmin(distances)
            central_points.append(cluster_points[central_point_idx])

        # Centroids GeoDataFrame
        central_features = [
            {
                "geometry": Point(central_point),
                "cluster": i,
                "name_microgrid": grid_name,
            }
            for i, central_point in enumerate(central_points)
        ]
        central_features_gdf = gpd.GeoDataFrame(
            central_features, crs=filtered_buildings.crs
        ).to_crs(geo_crs)
        all_central_features = pd.concat(
            [all_central_features, central_features_gdf], ignore_index=True
        )

        # Assign clusters to buildings
        filtered_buildings = filtered_buildings.copy()
        filtered_buildings["cluster_id"] = kmeans.labels_.astype(int)
        all_microgrid_buildings = pd.concat(
            [all_microgrid_buildings, filtered_buildings], ignore_index=True
        )

        # Counts per cluster
        buildings_class = (
            filtered_buildings.groupby("cluster_id")
            .tags_building.value_counts()
            .reset_index(name="count")
        )
        buildings_class["name_microgrid"] = grid_name
        all_buildings_class = pd.concat(
            [all_buildings_class, buildings_class], ignore_index=True
        )

    # Save centroids (unchanged)
    all_central_features.to_file(output_filepath_centroids, driver="GeoJSON")

    # Ensure unique columns, save buildings with GeoPandas (no manual JSON)
    if all_microgrid_buildings.columns.duplicated().any():
        all_microgrid_buildings = all_microgrid_buildings.loc[
            :, ~all_microgrid_buildings.columns.duplicated(keep="first")
        ].copy()

    gdf_out = all_microgrid_buildings.to_crs(geo_crs)
    gdf_out.to_file(output_filepath_buildings, driver="GeoJSON")

    # Save CSV (unchanged)
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
    geo_crs = snakemake.params.crs["geo_crs"]

    get_central_points_geojson_with_buildings(
        snakemake.input["buildings_geojson"],
        snakemake.output["clusters"],
        snakemake.config["buildings"]["n_clusters"],
        crs,
        house_area_limit,
        snakemake.output["clusters_with_buildings"],
        snakemake.output["buildings_type"],
        snakemake.config["microgrids_list"],
        geo_crs,
    )

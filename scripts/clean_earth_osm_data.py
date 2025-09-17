# -*- coding: utf-8 -*-

# TODO: Add docstring

import json
import os
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _helpers_dist import configure_logging, sets_path_to_root
from shapely.geometry import Point, Polygon


def extract_points(microgrid_shape_path, buildings_path, output_path):
    """
    From the downloaded data, extracts buildings located within the boundaries of each microgrid geometry
    and associates them with the respective microgrid name.

    Parameters
    ----------
    microgrid_shape_path : str
        Path to the GeoJSON file containing microgrid geometries.
    buildings_path : str
        Path to the GeoJSON file containing building geometries.
    output_path : str
        Path where the resulting GeoJSON file will be saved.

    Returns
    -------
    GeoDataFrame
        A GeoDataFrame containing the filtered buildings with an added field "name_microgrid"
        that associates each building to its corresponding microgrid.
    """
    # Load the GeoJSON files
    microgrid = gpd.read_file(microgrid_shape_path)
    buildings = gpd.read_file(buildings_path)
    # Create a GeoDataFrame to accumulate the results
    result = gpd.GeoDataFrame(columns=buildings.columns)
    # Iterate over each microgrid geometry
    for idx, microgrid_shape in microgrid.iterrows():
        # Extract the name of the microgrid
        microgrid_name = microgrid_shape["name"]
        # Filter buildings located within the microgrid geometry
        buildings_in_microgrid = buildings[
            buildings.geometry.within(microgrid_shape.geometry)
        ]
        # Add or replace the "name_microgrid" field with the microgrid name
        buildings_in_microgrid = buildings_in_microgrid.copy()
        buildings_in_microgrid["name_microgrid"] = microgrid_name
        # Append the filtered buildings to the final result
        result = gpd.GeoDataFrame(
            pd.concat([result, buildings_in_microgrid], ignore_index=True)
        )
    geojson_features = []
    try:
        for feat in result.iterfeatures(na="drop", drop_id=True):
            geojson_features.append(
                json.dumps(feat, ensure_ascii=False, separators=(",", ":"))
            )
    except TypeError:
        raw = json.loads(result.to_json())
        for feat in raw.get("features", []):
            props = feat.get("properties", {}) or {}
            for k in [k for k, v in list(props.items()) if v is None]:
                props.pop(k, None)
            geojson_features.append(
                json.dumps(feat, ensure_ascii=False, separators=(",", ":"))
            )
    outpath = Path(output_path)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "w", encoding="utf-8") as f:
        f.write('{"type":"FeatureCollection","features":[\n')
        f.write(",\n".join(geojson_features))
        f.write("\n]}\n")
    return result


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

# -*- coding: utf-8 -*-

# TODO: Add docstring

import json
import os

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _helpers_dist import configure_logging, sets_path_to_root
from shapely.geometry import Point, Polygon


def extract_points(microgrid_shape_path, buildings_path, output_path):
    # Carica i file GeoJSON
    microgrid = gpd.read_file(microgrid_shape_path)
    buildings = gpd.read_file(buildings_path)

    # Crea un GeoDataFrame per accumulare i risultati
    result = gpd.GeoDataFrame(columns=buildings.columns)

    # Itera su ogni geometria della microrete
    for idx, microgrid_shape in microgrid.iterrows():
        # Estrai il nome della microrete
        microgrid_name = microgrid_shape["name"]

        # Filtra gli edifici che si trovano nella geometria della microrete
        buildings_in_microgrid = buildings[buildings.geometry.within(microgrid_shape.geometry)]

        # Aggiungi o sostituisci il campo "name_microgrid" con il nome calcolato
        buildings_in_microgrid = buildings_in_microgrid.copy()
        buildings_in_microgrid["name_microgrid"] = microgrid_name

        # Aggiungi gli edifici filtrati al risultato finale
        result = gpd.GeoDataFrame(pd.concat([result, buildings_in_microgrid], ignore_index=True))

    # Salva il risultato come GeoJSON
    result.to_file(output_path, driver="GeoJSON")

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

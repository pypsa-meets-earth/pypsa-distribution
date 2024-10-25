# -*- coding: utf-8 -*-
"""
Estimates the population and the electric load of each microgrid.

Relevant Settings
-----------------

.. code:: yaml

    microgrids_list:
        microgridX: 
          lon_min:
          lon_max: 
          lat_min: 
          lat_max: 
    load:
        scaling_factor:

Inputs
------
- ``data/sample_profile.csv``: a load profile, which will be scaled through a scaling_factor to obtain the per person load

Outputs
-------
- ``resources/shapes/microgrid_shapes.geojson``: a geojson file of the shape of each microgrid,
- ``resources/masked_files/masked_file_{i+1}.tif``,
- ``resources/demand/microgrid_load_{i+1}.csv``: the electric load of the microgid,

Description
-----------
The rule :mod:`build_demand` contains functions that are used to create a shape file of the microgrid, to mask a raster with the shape file and to estimate 
the population. Then the population is multiplied for the per person load and the microgrid load is then obtained. The process applies to all the microgrids specified in config.yaml.
"""

import json
import logging
import os
import shutil

import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
import rasterio
import rasterio.mask
import requests
from _helpers_dist import (
    configure_logging,
    sets_path_to_root,
    two_2_three_digits_country,
)

_logger = logging.getLogger(__name__)
_logger.setLevel(logging.INFO)


def get_WorldPop_data(
    country_code,
    year,
    update=False,
    out_logging=False,
    size_min=300,
):
    """
    Download tiff file for each country code using the standard method from worldpop datastore with 1kmx1km resolution.

    Parameters
    ----------
    country_code : str
        Two letter country codes of the downloaded files.
        Files downloaded from https://data.worldpop.org/ datasets WorldPop UN adjusted
    year : int
        Year of the data to download
    update : bool
        Update = true, forces re-download of files
    size_min : int
        Minimum size of each file to download
    Returns
    -------
    WorldPop_inputfile : str
        Path of the file
    """

    three_digits_code = two_2_three_digits_country(country_code)

    if out_logging:
        _logger.info("Get WorldPop datasets")

    if country_code == "XK":
        WorldPop_filename = f"srb_ppp_{year}_UNadj_constrained.tif"
        WorldPop_urls = [
            f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/BSGM/SRB/{WorldPop_filename}",
            f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/maxar_v1/SRB/{WorldPop_filename}",
        ]
    else:
        WorldPop_filename = (
            f"{three_digits_code.lower()}_ppp_{year}_UNadj_constrained.tif"
        )
        # Urls used to possibly download the file
        WorldPop_urls = [
            f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/BSGM/{two_2_three_digits_country(country_code).upper()}/{WorldPop_filename}",
            f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/maxar_v1/{two_2_three_digits_country(country_code).upper()}/{WorldPop_filename}",
        ]

    WorldPop_inputfile = os.path.join(
        os.getcwd(),
        "pypsa-earth",
        "data",
        "WorldPop",
        WorldPop_filename,
    )  # Input filepath tif

    if not os.path.exists(WorldPop_inputfile) or update is True:
        if out_logging:
            _logger.warning(
                f"{WorldPop_filename} does not exist, downloading to {WorldPop_inputfile}"
            )
        #  create data/osm directory
        os.makedirs(os.path.dirname(WorldPop_inputfile), exist_ok=True)

        loaded = False
        for WorldPop_url in WorldPop_urls:
            with requests.get(WorldPop_url, stream=True) as r:
                with open(WorldPop_inputfile, "wb") as f:
                    if float(r.headers["Content-length"]) > size_min:
                        shutil.copyfileobj(r.raw, f)
                        loaded = True
                        break
        if not loaded:
            _logger.error(f"Impossible to download {WorldPop_filename}")

    return WorldPop_inputfile, WorldPop_filename


# Estimate the total population of tghe microgrid
def estimate_microgrid_population(
    n, p, raster_path, shapes_path, sample_profile, output_file
):
    # Read the sample profile of electricity demand and extract the column corresponding to the electric load
    per_unit_load = pd.read_csv(sample_profile)["0"] / p

    # Dataframe of the load
    microgrid_load = pd.DataFrame()

    # Load the GeoJSON file with the shapes to mask the raster
    shapes = gpd.read_file(shapes_path)

    # Mask the raster with each shape and save each masked raster as a new file
    for i, shape in shapes.iterrows():
        with rasterio.open(raster_path) as src:
            # Mask the raster with the current shape
            masked, out_transform = rasterio.mask.mask(src, [shape.geometry], crop=True)
            out_meta = src.meta.copy()
            out_meta.update(
                {
                    "driver": "GTiff",
                    "height": masked.shape[1],
                    "width": masked.shape[2],
                    "transform": out_transform,
                }
            )

        pop_microgrid = masked[masked >= 0].sum()

        col_name = "microgrid_1_bus_572666767"
        microgrid_load[col_name] = per_unit_load * pop_microgrid

    return pop_microgrid, microgrid_load


def calculate_load(
    n,
    p,
    raster_path,
    shapes_path,
    sample_profile,
    geojson_file,
    output_file,
    input_path,
):
    # Estimate the microgrid population and load using the existing function
    pop_microgrid, microgrid_load = estimate_microgrid_population(
        n, p, raster_path, shapes_path, sample_profile, output_file
    )
    building_class = pd.read_csv(input_path)
    total_buildings = building_class["count"].sum()
    building_for_cluster = pd.DataFrame(
        building_class.groupby("cluster_id").sum()["count"]
    )
    population_per_building = pop_microgrid / total_buildings
    population_per_cluster = building_for_cluster * population_per_building
    per_unit_load = pd.read_csv(sample_profile)["0"] / p
    load_per_cluster = population_per_cluster["count"].apply(
        lambda x: x * per_unit_load
    )
    load_per_cluster = load_per_cluster.T
    new_column_names = {i: f"bus_{i}" for i in range(load_per_cluster.shape[1])}
    load_per_cluster.rename(columns=new_column_names, inplace=True)
    load_per_cluster.insert(0, "snapshots", n.snapshots)
    load_per_cluster.to_csv(output_file, index=True)

    return load_per_cluster


def calculate_load_ramp(
    input_file_buildings,
    n,
    p,
    raster_path,
    shapes_path,
    sample_profile,
    output_file,
    input_file_profile_tier1,
    input_file_profile_tier2,
    input_file_profile_tier3,
    input_file_profile_tier4,
    input_file_profile_tier5,
    output_path_csv,
    tier_percent,
    date_start,
    date_end,
    inclusive,
):
    # Caricamento dei dati e calcolo della densità di popolazione
    cleaned_buildings = gpd.read_file(input_file_buildings)
    house = cleaned_buildings[cleaned_buildings["tags_building"] == "house"]
    pop_microgrid, microgrid_load = estimate_microgrid_population(
        n, p, raster_path, shapes_path, sample_profile, output_file
    )
    density = pop_microgrid / house["area_m2"].sum()

    # Calcolo superficie e popolazione per cluster
    grouped_buildings = cleaned_buildings.groupby("cluster_id")
    clusters = np.sort(cleaned_buildings["cluster_id"].unique())
    house_area_for_cluster = [
        grouped_buildings.get_group(cluster)[grouped_buildings.get_group(cluster)["tags_building"] == "house"]["area_m2"].sum()
        for cluster in clusters
    ]
    population_df = pd.DataFrame({"cluster": clusters, "house_area_for_cluster": house_area_for_cluster}).set_index("cluster")
    population_df["people_for_cluster"] = (population_df["house_area_for_cluster"] * density).round()
    tier_pop_df = pd.DataFrame(np.outer(population_df["people_for_cluster"], tier_percent), index=population_df.index)

    # Caricamento e creazione di DataFrames di domanda media e deviazione standard per ogni tier
    demand_files = [input_file_profile_tier1, input_file_profile_tier2, input_file_profile_tier3, input_file_profile_tier4, input_file_profile_tier5]
    mean_demand_tier_df = pd.DataFrame({f"tier_{i+1}": pd.read_excel(file)["mean"] for i, file in enumerate(demand_files)})
    mean_demand_tier_df.insert(0, "tier_0", np.zeros(len(mean_demand_tier_df)))
    mean_demand_tier_df.index = pd.date_range("00:00:00", periods=len(mean_demand_tier_df), freq="H").time

    if inclusive == "left":
        date_range = pd.date_range(start=date_start, end=date_end, freq="D")[:-1]
    else: 
        date_range = pd.date_range(start=date_start, end=date_end, freq="D")

    mean_demand_tier_df_extended = pd.concat([mean_demand_tier_df] * len(date_range), ignore_index=True)

    # Calcolo del carico totale per ogni cluster e tier
    result_dict = {}
    for k, pop_cluster in tier_pop_df.iterrows():
        load_df = pd.DataFrame()
        for j, n_person in enumerate(pop_cluster / 7):  # Scala la popolazione per famiglia
            mean_load = mean_demand_tier_df_extended.iloc[:, j] * n_person
            total_load = (mean_load) / 1e6
            load_df[f"tier_{j}"] = total_load
        result_dict[f"bus_{k}"] = load_df

    # Aggregazione del carico totale per cluster
    tot_result_dict = {f"{k}": df.sum(axis=1).rename(f"{k}") for k, df in result_dict.items()}
    tot_loads_df = pd.concat(tot_result_dict.values(), axis=1)
    if inclusive == "left":
        date_range_tot = pd.date_range(start=date_start, end=date_end, freq="H")[:-1]
    else: 
        date_range_tot = pd.date_range(start=date_start, end=date_end, freq="H")
    tot_loads_df.index=date_range_tot

    # Sostituzione dei valori zero con un valore minimo per evitare problemi di plotting
    small_value = 1e-26
    tot_loads_df.loc[:, (tot_loads_df == 0).all()] = small_value

    # Esportazione del DataFrame finale
    tot_loads_df.to_csv(output_path_csv)


def calculate_load_ramp_std(
    input_file_buildings,
    n,
    p,
    raster_path,
    shapes_path,
    sample_profile,
    output_file,
    input_file_profile_tier1,
    input_file_profile_tier2,
    input_file_profile_tier3,
    input_file_profile_tier4,
    input_file_profile_tier5,
    output_path_csv,
    tier_percent,
    date_start,
    date_end,
    inclusive,
):
    # Caricamento dei dati e calcolo della densità di popolazione
    cleaned_buildings = gpd.read_file(input_file_buildings)
    house = cleaned_buildings[cleaned_buildings["tags_building"] == "house"]
    pop_microgrid, microgrid_load = estimate_microgrid_population(
        n, p, raster_path, shapes_path, sample_profile, output_file
    )
    density = pop_microgrid / house["area_m2"].sum()

    # Calcolo superficie e popolazione per cluster
    grouped_buildings = cleaned_buildings.groupby("cluster_id")
    clusters = np.sort(cleaned_buildings["cluster_id"].unique())
    house_area_for_cluster = [
        grouped_buildings.get_group(cluster)[grouped_buildings.get_group(cluster)["tags_building"] == "house"]["area_m2"].sum()
        for cluster in clusters
    ]
    population_df = pd.DataFrame({"cluster": clusters, "house_area_for_cluster": house_area_for_cluster}).set_index("cluster")
    population_df["people_for_cluster"] = (population_df["house_area_for_cluster"] * density).round()
    tier_pop_df = pd.DataFrame(np.outer(population_df["people_for_cluster"], tier_percent), index=population_df.index)

    # Caricamento e creazione di DataFrames di domanda media e deviazione standard per ogni tier
    demand_files = [input_file_profile_tier1, input_file_profile_tier2, input_file_profile_tier3, input_file_profile_tier4, input_file_profile_tier5]
    mean_demand_tier_df = pd.DataFrame({f"tier_{i+1}": pd.read_excel(file)["mean"] for i, file in enumerate(demand_files)})
    std_demand_tier_df = pd.DataFrame({f"tier_{i+1}": pd.read_excel(file)["std"] for i, file in enumerate(demand_files)})
    mean_demand_tier_df.insert(0, "tier_0", np.zeros(len(mean_demand_tier_df)))
    std_demand_tier_df.insert(0, "tier_0", np.zeros(len(mean_demand_tier_df)))
    mean_demand_tier_df.index = pd.date_range("00:00:00", periods=len(mean_demand_tier_df), freq="H").time
    std_demand_tier_df.index = pd.date_range("00:00:00", periods=len(mean_demand_tier_df), freq="H").time

    if inclusive == "left":
        date_range = pd.date_range(start=date_start, end=date_end, freq="D")[:-1]
    else: 
        date_range = pd.date_range(start=date_start, end=date_end, freq="D")

    mean_demand_tier_df_extended = pd.concat([mean_demand_tier_df] * len(date_range), ignore_index=True)
    std_demand_tier_df_extended = pd.concat([std_demand_tier_df] * len(date_range), ignore_index=True)

    # Calcolo del carico totale per ogni cluster e tier
    result_dict = {}
    for k, pop_cluster in tier_pop_df.iterrows():
        load_df = pd.DataFrame()
        for j, n_person in enumerate(pop_cluster / 7):  # Scala la popolazione per famiglia
            mean_load = mean_demand_tier_df_extended.iloc[:, j] * n_person
            std_load = np.random.normal(mean_demand_tier_df_extended.iloc[:, j], std_demand_tier_df_extended.iloc[:, j]) * np.sqrt(n_person)
            total_load = (mean_load + std_load) / 1e6
            load_df[f"tier_{j}"] = total_load
        result_dict[f"bus_{k}"] = load_df

    # Aggregazione del carico totale per cluster
    tot_result_dict = {f"{k}": df.sum(axis=1).rename(f"{k}") for k, df in result_dict.items()}
    tot_loads_df = pd.concat(tot_result_dict.values(), axis=1)
    if inclusive == "left":
        date_range_tot = pd.date_range(start=date_start, end=date_end, freq="H")[:-1]
    else: 
        date_range_tot = pd.date_range(start=date_start, end=date_end, freq="H")
    tot_loads_df.index=date_range_tot

    # Sostituzione dei valori zero con un valore minimo per evitare problemi di plotting
    small_value = 1e-26
    tot_loads_df.loc[:, (tot_loads_df == 0).all()] = small_value

    # Esportazione del DataFrame finale
    tot_loads_df.to_csv(output_path_csv)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_dist import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("build_demand")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.create_network)
    sample_profile = snakemake.input["sample_profile"]
    tier_percent = snakemake.params.tier["tier_percent"]
    date_start=snakemake.params.snapshots["start"]
    date_end=snakemake.params.snapshots["end"]
    inclusive=snakemake.params.snapshots["inclusive"]

    build_demand_model = snakemake.params.build_demand_model["type"]

    assert (
        len(snakemake.config["countries"]) == 1
    ), "Error: only a country shall be specified"

    worldpop_path, worldpop_flname = get_WorldPop_data(
        snakemake.config["countries"][
            0
        ],  # TODO: this needs fix to generalize the countries
        snakemake.config["build_shape_options"]["year"],
        False,
    )

    estimate_microgrid_population(
        n,
        snakemake.config["load"]["scaling_factor"],
        worldpop_path,
        snakemake.input["microgrid_shapes"],
        sample_profile,
        snakemake.output["electric_load"],
    )
    if build_demand_model == 0:
        calculate_load(
            n,
            snakemake.config["load"]["scaling_factor"],
            worldpop_path,
            snakemake.input["microgrid_shapes"],
            sample_profile,
            snakemake.input["clusters_with_buildings"],
            snakemake.output["electric_load"],
            snakemake.input["building_csv"],
        )

    elif build_demand_model == 1:
        calculate_load_ramp(
            snakemake.input["clusters_with_buildings"],
            n,
            snakemake.config["load"]["scaling_factor"],
            worldpop_path,
            snakemake.input["microgrid_shapes"],
            sample_profile,
            snakemake.output["electric_load"],
            snakemake.input["profile_Tier1"],
            snakemake.input["profile_Tier2"],
            snakemake.input["profile_Tier3"],
            snakemake.input["profile_Tier4"],
            snakemake.input["profile_Tier5"],
            snakemake.output["electric_load"],
            tier_percent,
            date_start,
            date_end,
            inclusive,
        )
    elif build_demand_model == 2:

        calculate_load_ramp_std(
            snakemake.input["clusters_with_buildings"],
            n,
            snakemake.config["load"]["scaling_factor"],
            worldpop_path,
            snakemake.input["microgrid_shapes"],
            sample_profile,
            snakemake.output["electric_load"],
            snakemake.input["profile_Tier1"],
            snakemake.input["profile_Tier2"],
            snakemake.input["profile_Tier3"],
            snakemake.input["profile_Tier4"],
            snakemake.input["profile_Tier5"],
            snakemake.output["electric_load"],
            tier_percent,
            date_start,
            date_end,
            inclusive,
        )

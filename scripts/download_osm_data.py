# -*- coding: utf-8 -*-
import logging
import os
import shutil
from pathlib import Path

import yaml
from _helpers_dist import configure_logging, create_logger, read_osm_config
from earth_osm import eo

logger = create_logger(__name__)


def country_list_to_geofk(country_list):
    """
    Convert the requested country list into geofk norm.

    Parameters
    ----------
    input : str
        Any two-letter country name or aggregation of countries given in osm_config.yaml
        Country name duplications won't distort the result.
        Examples are:
        ["NG","ZA"], downloading osm data for Nigeria and South Africa
        ["SNGM"], downloading data for Senegal&Gambia shape
        ["NG","ZA","NG"], won't distort result.

    Returns
    -------
    full_codes_list : list
        Example ["NG","ZA"]
    """
    full_codes_list = [convert_iso_to_geofk(c_code) for c_code in set(country_list)]

    return full_codes_list


def convert_iso_to_geofk(
    iso_code, iso_coding=True, convert_dict=read_osm_config("iso_to_geofk_dict")
):
    """
    Function to convert the iso code name of a country into the corresponding
    geofabrik In Geofabrik, some countries are aggregated, thus if a single
    country is requested, then all the agglomeration shall be downloaded For
    example, Senegal (SN) and Gambia (GM) cannot be found alone in geofabrik,
    but they can be downloaded as a whole SNGM.

    The conversion directory, initialized to iso_to_geofk_dict is used to perform such conversion
    When a two-letter code country is found in convert_dict, and iso_coding is enabled,
    then that two-letter code is converted into the corresponding value of the dictionary

    Parameters
    ----------
    iso_code : str
        Two-code country code to be converted
    iso_coding : bool
        When true, the iso to geofk is performed
    convert_dict : dict
        Dictionary used to apply the conversion iso to geofk
        The keys correspond to the countries iso codes that need a different region to be downloaded
    """
    if iso_coding and iso_code in convert_dict:
        return convert_dict[iso_code]
    else:
        return iso_code


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_dist import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("download_osm_data")
        sets_path_to_root("pypsa-distribution")
    configure_logging(snakemake)

    run = snakemake.config.get("run", {})
    RDIR = run["name"] + "/" if run.get("name") else ""
    store_path_resources = Path.cwd() / "resources" / RDIR / "osm" / "raw"
    store_path_data = Path.cwd() / "data" / "osm"
    countries = snakemake.config["countries"]
    country_list = country_list_to_geofk(countries)

    eo.save_osm_data(
        region_list=country_list,
        primary_name="building",
        feature_list=["ALL"],
        update=False,
        mp=False,
        data_dir=store_path_data,
        out_dir=store_path_resources,
        out_format=["csv", "geojson"],
        out_aggregate=True,
    )

    out_path = Path.joinpath(store_path_resources, "out")
    out_formats = ["csv", "geojson"]
    new_files = os.listdir(out_path)

    for f in out_formats:
        new_file_name = Path.joinpath(store_path_resources, f"all_raw_building.{f}")
        old_file = list(Path(out_path).glob(f"*building.{f}"))

        if not old_file:
            with open(new_file_name, "w") as f:
                pass
        else:
            logger.info(f"Move {old_file[0]} to {new_file_name}")
            shutil.move(old_file[0], new_file_name)

# -*- coding: utf-8 -*-
import json
import logging
import os
import shutil
from pathlib import Path
import geopandas as gpd
import mercantile
import pandas as pd
import requests
from _helpers_dist import configure_logging, create_logger, read_osm_config
from earth_osm import eo
from shapely import geometry
from shapely.geometry import mapping
from tqdm import tqdm

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


def retrieve_osm_data_geojson(microgrids_list, features, url, path):
    """
    Downloads OSM data for buildings, low/medium voltage lines, and distribution network components
    from the Overpass API, for each microgrid defined in the configuration.

    The function queries Overpass using bounding boxes for each microgrid and extracts relevant
    OSM features. These features are converted into GeoJSON format and saved as separate .geojson
    files by feature type (e.g., buildings, power lines, substations).

    Parameters
    ----------
    microgrids_list : dict
        Dictionary of microgrid definitions with bounding box coordinates present in the config. file.

    features : list of str
        List of OSM feature types to extract. Supported values:
        - "building"
        - "minor_line"
        - "generator"
        - "substation_and_pole"

    url : str
        Overpass API endpoint URL

    path : str
        Directory path where the resulting .geojson files will be saved.
    """

    for feature in features:
        geojson_features = []

        for grid_name, grid_data in microgrids_list.items():
            lat_min = grid_data["lat_min"]
            lon_min = grid_data["lon_min"]
            lat_max = grid_data["lat_max"]
            lon_max = grid_data["lon_max"]

            # Determina filename e geometry_type UNA VOLTA per ogni feature_name
            if feature == "building":
                filename = "all_raw_buildings.geojson"
                geometry_type = "Polygon"
                overpass_query = f"""
                [out:json][timeout:60];
                (
                way["building"]({lat_min},{lon_min},{lat_max},{lon_max});
                );
                (._;>;);
                out body;
                """

            elif feature == "minor_line":
                filename = "all_raw_line.geojson"
                geometry_type = "LineString"
                overpass_query = f"""
                [out:json][timeout:60];
                (
                way["power"~"^(cable|minor_line)$"]({lat_min},{lon_min},{lat_max},{lon_max});
                relation["power"~"^(cable|minor_line)$"]({lat_min},{lon_min},{lat_max},{lon_max});
                );
                (._;>;);
                out body;
                """

            elif feature == "generator":
                filename = "all_raw_generator.geojson"
                geometry_type = "Polygon"
                overpass_query = f"""
                [out:json][timeout:60];
                (
                node["power"~"generator|plant"]({lat_min},{lon_min},{lat_max},{lon_max});
                way["power"~"generator|plant"]({lat_min},{lon_min},{lat_max},{lon_max});
                relation["power"~"generator|plant"]({lat_min},{lon_min},{lat_max},{lon_max});
                );
                (._;>;);
                out body;
                """

            elif feature == "substation_and_pole":
                filename = "all_raw_substation.geojson"
                geometry_type = "Polygon"
                overpass_query = f"""
                [out:json][timeout:60];
                (
                node["power"="pole"]({lat_min},{lon_min},{lat_max},{lon_max});
                node["power"="substation"]({lat_min},{lon_min},{lat_max},{lon_max});
                way["power"="substation"]({lat_min},{lon_min},{lat_max},{lon_max});
                relation["power"="substation"]({lat_min},{lon_min},{lat_max},{lon_max});
                );
                (._;>;);
                out body;
                """

            else:
                logger.error(f"Unsupported feature: {feature}")
                continue

            try:
                logger.info(
                    f"Querying Overpass API for microgrid: {grid_name} with feature: {feature}"
                )
                response = requests.get(url, params={"data": overpass_query})
                response.raise_for_status()
                data = response.json()

                if "elements" not in data:
                    logger.error(
                        f"No elements found for microgrid: {grid_name} with feature: {feature}"
                    )
                    continue

                node_coordinates = {
                    node["id"]: [node["lon"], node["lat"]]
                    for node in data["elements"]
                    if node["type"] == "node"
                }

                way_elements = {
                    element["id"]: element
                    for element in data["elements"]
                    if element["type"] == "way" and "nodes" in element
                }

                for element in data["elements"]:
                    if element["type"] == "way" and "nodes" in element:
                        logger.debug(f"Processing {feature}: {element['id']}")
                        coordinates = [
                            node_coordinates[node_id]
                            for node_id in element["nodes"]
                            if node_id in node_coordinates
                        ]

                        if not coordinates:
                            logger.warning(
                                f"No coordinates for {feature}: {element['id']}"
                            )
                            continue

                        properties = {"name_microgrid": grid_name, "id": element["id"]}
                        if "tags" in element:
                            properties.update(element["tags"])

                        feature = {
                            "type": "Feature",
                            "properties": properties,
                            "geometry": {
                                "type": geometry_type,
                                "coordinates": (
                                    [coordinates]
                                    if geometry_type == "Polygon"
                                    else coordinates
                                ),
                            },
                        }
                        geojson_features.append(
                            json.dumps(feature, separators=(",", ":"))
                        )

            except json.JSONDecodeError:
                logger.error(
                    f"JSON decoding error for microgrid: {grid_name} with feature: {feature}"
                )
            except requests.exceptions.RequestException as e:
                logger.error(
                    f"Request error for microgrid: {grid_name}: {e} with feature: {feature}"
                )

        try:
            outpath = Path(path) / filename
            outpath.parent.mkdir(parents=True, exist_ok=True)

            with open(outpath, "w") as f:
                f.write('{"type":"FeatureCollection","features":[\n')
                f.write(
                    ",\n".join(geojson_features)
                )  # Write features in one-line format
                f.write("\n]}\n")

            logger.info(f"GeoJSON saved to {outpath}")

        except IOError as e:
            logger.error(f"Error saving GeoJSON file: {e}")


def retrive_and_merge_osm_with_ml(microgrid_list, url, osm_path, export_path):

    link = pd.read_csv(url, dtype=str)
    mML_gdf = gpd.GeoDataFrame()
    idx = 0

    for gridname, grid_data in microgrid_list.items():
        lon_min, lat_min, lon_max, lat_max = (
            grid_data["lon_min"],
            grid_data["lat_min"],
            grid_data["lon_max"],
            grid_data["lat_max"],
        )
        microgrid_shape = geometry.box(lon_min, lat_min, lon_max, lat_max)

        quad_keys = list({
            mercantile.quadkey(tile)
            for tile in mercantile.tiles(lon_min, lat_min, lon_max, lat_max, zooms=9)
        })

        logger.info(f"[{gridname}] AOI spans {len(quad_keys)} tiles.")

        for quad_key in quad_keys:
            row = link[link["QuadKey"] == quad_key]

            if row.shape[0] == 1:
                json_url = row.iloc[0]["Url"]
                logger.info(f"[{gridname}] Downloading {json_url} for quad_key {quad_key}")
                try:
                    df_json = pd.read_json(json_url, lines=True)
                except Exception as e:
                    logger.warning(f"[{gridname}] Failed to download {json_url}: {e}")
                    continue

                if "geometry" not in df_json.columns:
                    logger.warning(f"[{gridname}] No geometry found for quad_key {quad_key}")
                    continue

                df_json["geometry"] = df_json["geometry"].apply(geometry.shape)
                gdf = gpd.GeoDataFrame(df_json, geometry="geometry", crs="EPSG:4326")
                gdf = gdf[gdf.geometry.within(microgrid_shape)]
                gdf["id"] = range(idx, idx + len(gdf))
                idx += len(gdf)
                mML_gdf = pd.concat([mML_gdf, gdf], ignore_index=True)

            elif row.shape[0] > 1:
                logger.error(f"Multiple entries found for QuadKey: {quad_key}")
                continue
            else:
                logger.info(f"QuadKey not found: {quad_key}")
                continue

    if mML_gdf.empty:
        logger.warning("No Microsoft ML buildings found in the AOI.")
        return

    mML_gdf["height"] = mML_gdf["properties"].apply(
        lambda s: s.get("height") if isinstance(s, dict) else None
    )

    osm_gdf = gpd.read_file(osm_path)
    if "geometry" not in osm_gdf.columns:
        raise ValueError("OSM file does not contain 'geometry' column.")

    crs_meters = "EPSG:32632"
    mML_gdf = mML_gdf.to_crs(crs_meters)
    osm_gdf = osm_gdf.to_crs(crs_meters)

    mML_buffered = mML_gdf.copy()
    mML_buffered["geometry"] = mML_buffered.geometry.buffer(1.0)

    joined = gpd.sjoin(osm_gdf, mML_buffered, how="inner", predicate="intersects")
    joined = joined.astype({"index_right": int})
    joined = joined.drop(columns=["properties"], errors="ignore")

    merged_features = []
    for _, row in joined.iterrows():
        idx = row["index_right"]
        mML_row = mML_gdf.loc[idx]

        osm_props = {
            k: v
            for k, v in row.drop(labels=["index_right", "geometry"]).items()
            if pd.notna(v)
        }

        if "height" in mML_row and pd.notna(mML_row["height"]):
            osm_props["height"] = mML_row["height"]

        geometry_ = row.geometry
        merged_features.append({"geometry": geometry_, **osm_props})

    merged_duplicates_gdf = gpd.GeoDataFrame(merged_features, crs=crs_meters)
    duplicate_idxs = joined.index
    cluster_unique = osm_gdf.drop(index=duplicate_idxs)
    final_gdf = pd.concat([merged_duplicates_gdf, cluster_unique], ignore_index=True)
    final_gdf = final_gdf.to_crs("EPSG:4326")

    features = []

    for _, row in final_gdf.iterrows():
        props = {
            k: v
            for k, v in row.items()
            if k != "geometry"
            and pd.notna(v)
            and str(v).strip().lower() not in ["", "null", "none"]
        }
        geom = row["geometry"]
        if geom is not None and not geom.is_empty:
            features.append({
                "type": "Feature",
                "properties": props,
                "geometry": mapping(geom)
            })

    with open(export_path, "w") as f:
        f.write('{"type":"FeatureCollection","features":[\n')
        for i, feature in enumerate(features):
            line = json.dumps(feature, separators=(",", ":"))
            if i < len(features) - 1:
                f.write(line + ",\n")
            else:
                f.write(line + "\n")
        f.write("]}")

    logger.info(" Merge completato e salvato come GeoJSON formattato riga-per-riga.")




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

    if snakemake.config["enable"]["download_osm_method"] == "earth_osm":
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

    elif snakemake.config["enable"]["download_osm_method"] == "overpass":
        microgrids_list = snakemake.config["microgrids_list"]
        features = ["building", "minor_line", "generator", "substation_and_pole"]
        overpass_url = "https://overpass-api.de/api/interpreter"
        output_file = Path.cwd() / "resources" / RDIR / "osm" / "raw"
        retrieve_osm_data_geojson(microgrids_list, features, overpass_url, output_file)
        if snakemake.config["enable"]["download_and_merge_microsoft_ML_building"]:
            osm_path = (
                Path.cwd()
                / "resources"
                / RDIR
                / "osm"
                / "raw"
                / "all_raw_buildings.geojson"
            )
            microsoft_data_url = "https://minedbuildings.z5.web.core.windows.net/global-buildings/dataset-links.csv"
            export_path = (
                Path.cwd()
                / "resources"
                / RDIR
                / "osm"
                / "raw"
                / "all_raw_buildings.geojson"
            )
            retrive_and_merge_osm_with_ml(
                microgrids_list, microsoft_data_url, osm_path, export_path
            )

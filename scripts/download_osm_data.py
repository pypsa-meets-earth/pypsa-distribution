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
import yaml
from _helpers_dist import configure_logging, create_logger, read_osm_config
from earth_osm import eo
from shapely import geometry
from shapely.geometry import mapping

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


def retrieve_osm_data_geojson(microgrids_list, feature, url, path):
    """
    The buildings inside the specified coordinates are retrieved by using overpass API.
    The region coordinates should be defined in the config.yaml file.
    Parameters
    ----------
    microgrids_list : dict
        Dictionary containing the microgrid names and their bounding box coordinates (lat_min, lon_min, lat_max, lon_max).
    features : str
        The feature that is searched in the osm database
    url : str
        osm query address
    path : str
        Directory where the GeoJSON file will be saved.
    """
    # Collect all features from all microgrids
    for feature in features:
        geojson_features = []

        for grid_name, grid_data in microgrids_list.items():
            # Extract the bounding box coordinates for the current microgrid to construct the query
            lat_min = grid_data["lat_min"]
            lon_min = grid_data["lon_min"]
            lat_max = grid_data["lat_max"]
            lon_max = grid_data["lon_max"]

            if feature == "building":
                filename = "all_raw_building.geojson"
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
                way["power"~"^(minor_line)$"]({lat_min},{lon_min},{lat_max},{lon_max});
                relation["power"~"^(minor_line)$"]({lat_min},{lon_min},{lat_max},{lon_max});
                );
                (._;>;);
                out body;
                """
            
            elif feature == "cable":
                filename = "all_raw_cable.geojson"
                geometry_type = "LineString"
                overpass_query = f"""
                [out:json][timeout:60];
                (
                way["power"~"^(cable)$"]({lat_min},{lon_min},{lat_max},{lon_max});
                relation["power"~"^(cable)$"]({lat_min},{lon_min},{lat_max},{lon_max});
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
                )  # Log the current query
                response = requests.get(
                    url, params={"data": overpass_query}
                )  # Send the query to Overpass API
                response.raise_for_status()  # Raise an error if the request fails
                data = response.json()  # Parse the JSON response

                # Check if the response contains any elements
                if "elements" not in data:
                    logger.error(
                        f"No elements found for microgrid: {grid_name} with feature: {feature}"
                    )
                    continue
                # Extract node coordinates from the response
                node_coordinates = {
                    node["id"]: [node["lon"], node["lat"]]
                    for node in data["elements"]
                    if node["type"] == "node"
                }

                if feature == "substation_and_pole":
                    for element in data["elements"]:
                        if (
                            element["type"] == "node"
                            and "lon" in element
                            and "lat" in element
                        ):
                            properties = {
                                "name_microgrid": grid_name,
                                "id": element["id"],
                            }
                            tags = element.get("tags")
                            if isinstance(tags, dict):
                                properties.update(tags)

                            geometry = {
                                "type": "Point",
                                "coordinates": [element["lon"], element["lat"]],
                            }

                            feature = {
                                "type": "Feature",
                                "properties": properties,
                                "geometry": geometry,
                            }

                            geojson_features.append(
                                json.dumps(feature, separators=(",", ":"))
                            )
                        # Process "way" elements to construct polygon geometries
                        elif element["type"] == "way" and "nodes" in element:
                            # Get the coordinates of the nodes that form the way
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

                            # Add properties for the feature, including the microgrid name and element ID
                            properties = {
                                "name_microgrid": grid_name,
                                "id": element["id"],
                            }
                            tags = element.get("tags")
                            if isinstance(tags, dict):
                                properties.update(tags)

                            # Create a GeoJSON feature for the way
                            if geometry_type == "Polygon":
                                geometry = {
                                    "type": "Polygon",
                                    "coordinates": [coordinates],  # Close the polygon
                                }
                            else:
                                geometry = {
                                    "type": "LineString",
                                    "coordinates": coordinates,
                                }
                            feature = {
                                "type": "Feature",
                                "properties": properties,
                                "geometry": geometry,
                            }
                            # Serialize each feature as a compact JSON string and add it to the list
                            geojson_features.append(
                                json.dumps(feature, separators=(",", ":"))
                            )
                else:
                    for element in data["elements"]:
                        if element["type"] == "way" and "nodes" in element:
                            # Get the coordinates of the nodes that form the way
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

                            # Add properties for the feature, including the microgrid name and element ID
                            properties = {
                                "name_microgrid": grid_name,
                                "id": element["id"],
                            }
                            tags = element.get("tags")
                            if isinstance(tags, dict):
                                properties.update(tags)

                            # Create a GeoJSON feature for the way
                            if geometry_type == "Polygon":
                                geometry = {
                                    "type": "Polygon",
                                    "coordinates": [coordinates],  # Close the polygon
                                }
                            else:
                                geometry = {
                                    "type": "LineString",
                                    "coordinates": coordinates,
                                }
                            feature = {
                                "type": "Feature",
                                "properties": properties,
                                "geometry": geometry,
                            }
                            # Serialize each feature as a compact JSON string and add it to the list
                            geojson_features.append(
                                json.dumps(feature, separators=(",", ":"))
                            )

            except json.JSONDecodeError:
                # Handle JSON parsing errors
                logger.error(
                    f"JSON decoding error for microgrid: {grid_name}  with feature: {feature}"
                )
            except requests.exceptions.RequestException as e:
                # Handle request-related errors
                logger.error(
                    f"Request error for microgrid: {grid_name}: {e}  with feature: {feature}"
                )

                # Save all features to a single GeoJSON file
            try:
                outpath = Path(path) / filename
                outpath.parent.mkdir(parents=True, exist_ok=True)

                with open(outpath, "w") as f:
                    f.write('{"type":"FeatureCollection","features":[\n')
                    f.write(
                        ",\n".join(geojson_features)
                    )  # Write features in one-line format
                    f.write("\n]}\n")

                logger.info(f"Combined GeoJSON saved to {outpath}")

            except IOError as e:
                logger.error(f"Error saving GeoJSON file: {e}")


def download_and_merge_Microsoft_buildings(url, microgrid_list, osm_path, output_path):
    # Load tile-to-URL mapping from Microsoft building dataset
    link = pd.read_csv(url, dtype=str)
    mML_gdf = gpd.GeoDataFrame()
    idx = 0

    # Process each microgrid bounding box
    for gridname, grid_data in microgrid_list.items():
        lon_min, lat_min, lon_max, lat_max = (
            grid_data["lon_min"],
            grid_data["lat_min"],
            grid_data["lon_max"],
            grid_data["lat_max"],
        )
        microgrid_shape = geometry.box(lon_min, lat_min, lon_max, lat_max)

        # Get quadkeys covering the bounding box at zoom level 9
        quad_keys = list(
            {
                mercantile.quadkey(tile)
                for tile in mercantile.tiles(
                    lon_min, lat_min, lon_max, lat_max, zooms=9
                )
            }
        )

        # Download and filter building geometries per quadkey
        for quad_key in quad_keys:
            row = link[link["QuadKey"] == quad_key]

            if row.shape[0] == 1:
                json_url = row.iloc[0]["Url"]
                try:
                    df_json = pd.read_json(json_url, lines=True)
                except Exception:
                    continue

                if "geometry" not in df_json.columns:
                    continue

                # Convert JSON geometry to shapely shape and keep only those inside the microgrid
                df_json["geometry"] = df_json["geometry"].apply(geometry.shape)
                gdf = gpd.GeoDataFrame(df_json, geometry="geometry", crs="EPSG:4326")
                gdf = gdf[gdf.geometry.within(microgrid_shape)]
                gdf["id"] = range(idx, idx + len(gdf))
                idx += len(gdf)
                mML_gdf = pd.concat([mML_gdf, gdf], ignore_index=True)

    # Exit if no buildings found
    if mML_gdf.empty:
        print("vuoto")
        return

    # Extract height field from nested 'properties' dictionary
    mML_gdf["height"] = mML_gdf["properties"].apply(
        lambda s: (
            max(s.get("height", 0), 0)
            if isinstance(s, dict) and isinstance(s.get("height", None), (int, float))
            else 0
        )
    )

    # Load OSM buildings
    osm_gdf = gpd.read_file(osm_path)
    if "geometry" not in osm_gdf.columns:
        raise ValueError("OSM file does not contain 'geometry' column.")

    # Convert both datasets to a common metric CRS
    crs_meters = "EPSG:32632"
    mML_gdf = mML_gdf.to_crs(crs_meters)
    osm_gdf = osm_gdf.to_crs(crs_meters)

    # Slightly buffer Microsoft buildings for spatial matching
    mML_buffered = mML_gdf.copy()
    mML_buffered["geometry"] = mML_buffered.geometry.buffer(1.0)

    # Perform spatial join between OSM and buffered ML buildings
    joined = gpd.sjoin(osm_gdf, mML_buffered, how="inner", predicate="intersects")
    joined = joined.astype({"index_right": int})

    # Clean and merge matched attributes
    clean_join_rows = []
    for _, row in joined.iterrows():
        geometry_ = row["geometry"]
        cleaned_props = {
            k: v
            for k, v in row.items()
            if pd.notna(v) and v is not None and k not in ["geometry", "properties"]
        }
        clean_join_rows.append({**cleaned_props, "geometry": geometry_})

    cleaned_join_gdf = gpd.GeoDataFrame(clean_join_rows, crs=crs_meters)

    # Remove duplicates from original OSM set
    duplicate_idxs = joined.index
    cluster_unique = osm_gdf.drop(index=duplicate_idxs).to_crs(crs=crs_meters)

    # Merge cleaned duplicates with unique untouched buildings
    final_gdf = pd.concat([cleaned_join_gdf, cluster_unique], ignore_index=True)
    final_gdf = final_gdf.to_crs("EPSG:4326")

    # Convert merged GeoDataFrame to GeoJSON features
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
        safe_props = {
            k: (v.isoformat() if isinstance(v, pd.Timestamp) else v)
            for k, v in props.items()
        }
        if geom is not None and not geom.is_empty:
            features.append(
                {"type": "Feature", "properties": safe_props, "geometry": mapping(geom)}
            )

    # Write to GeoJSON file
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write('{"type":"FeatureCollection","features":[\n')
        f.write(",\n".join([json.dumps(feature) for feature in features]))
        f.write("\n]}\n")


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
        features = ["building", "minor_line", "cable", "generator", "substation_and_pole"]
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
                / "all_raw_building.geojson"
            )
            microsoft_data_url = "https://minedbuildings.z5.web.core.windows.net/global-buildings/dataset-links.csv"
            export_path = (
                Path.cwd()
                / "resources"
                / RDIR
                / "osm"
                / "raw"
                / "all_raw_building.geojson"
            )

            download_and_merge_Microsoft_buildings(
                microsoft_data_url, microgrids_list, osm_path, export_path
            )

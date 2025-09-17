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
    """
    full_codes_list = [convert_iso_to_geofk(c_code) for c_code in set(country_list)]
    return full_codes_list


def convert_iso_to_geofk(
    iso_code, iso_coding=True, convert_dict=read_osm_config("iso_to_geofk_dict")
):
    """
    Convert ISO country code to the geofabrik region code when needed.
    """
    if iso_coding and iso_code in convert_dict:
        return convert_dict[iso_code]
    else:
        return iso_code


import time
def retrieve_osm_data_geojson(microgrids_list, features, url, path):
    """
    Retrieve OpenStreetMap data from the Overpass API for a list of microgrids.
    For each requested OSM feature (buildings, lines, cables, generators,
    substations, poles), the function:
    • runs a single Overpass query per feature per microgrid bounding box,
    • converts the returned nodes/ways into GeoJSON Features with the proper
        geometry type (Polygon, LineString or Point),
    • and writes one GeoJSON file per feature into `path`.
    """

    for feature in features:
        geojson_features = []

        if feature == "building":
            filename = "all_raw_buildings.geojson"
            geometry_type = "Polygon"
        elif feature == "minor_line":
            filename = "all_raw_lines.geojson"
            geometry_type = "LineString"
        elif feature == "cable":
            filename = "all_raw_cables.geojson"
            geometry_type = "LineString"
        elif feature == "generator":
            filename = "all_raw_generators.geojson"
            geometry_type = "Polygon"
        elif feature == "substation":
            filename = "all_raw_substations.geojson"
            geometry_type = "Polygon"
        elif feature == "pole":
            filename = "all_raw_poles.geojson"
            geometry_type = "Point"
        else:
            logger.error(f"Unsupported feature: {feature}")
            continue

        for grid_name, grid_data in microgrids_list.items():
            lat_min = grid_data["lat_min"]; lon_min = grid_data["lon_min"]
            lat_max = grid_data["lat_max"]; lon_max = grid_data["lon_max"]

            
            if feature == "building":
                overpass_query = f"""
                [out:json][timeout:120];
                ( way["building"]({lat_min},{lon_min},{lat_max},{lon_max}); );
                (._;>;);
                out body;
                """
            elif feature == "minor_line":
                overpass_query = f"""
                [out:json][timeout:120];
                (
                  way["power"~"^(minor_line)$"]({lat_min},{lon_min},{lat_max},{lon_max});
                  relation["power"~"^(minor_line)$"]({lat_min},{lon_min},{lat_max},{lon_max});
                );
                (._;>;);
                out body;
                """
            elif feature == "cable":
                overpass_query = f"""
                [out:json][timeout:120];
                (
                  way["power"~"^(cable)$"]({lat_min},{lon_min},{lat_max},{lon_max});
                  relation["power"~"^(cable)$"]({lat_min},{lon_min},{lat_max},{lon_max});
                );
                (._;>;);
                out body;
                """
            elif feature == "generator":
                overpass_query = f"""
                [out:json][timeout:120];
                (
                  node["power"~"generator|plant"]({lat_min},{lon_min},{lat_max},{lon_max});
                  way["power"~"generator|plant"]({lat_min},{lon_min},{lat_max},{lon_max});
                  relation["power"~"generator|plant"]({lat_min},{lon_min},{lat_max},{lon_max});
                );
                (._;>;);
                out body;
                """
            elif feature == "substation":
                overpass_query = f"""
                [out:json][timeout:120];
                (
                  node["power"="substation"]({lat_min},{lon_min},{lat_max},{lon_max});
                  way["power"="substation"]({lat_min},{lon_min},{lat_max},{lon_max});
                );
                (._;>;);
                out body;
                """
            elif feature == "pole":
                overpass_query = f"""
                [out:json][timeout:120];
                (
                  node["power"="pole"]({lat_min},{lon_min},{lat_max},{lon_max});
                );
                (._;>;);
                out body;
                """

            try:
                logger.info(f"Querying Overpass API for microgrid: {grid_name} with feature: {feature}")
                response = requests.get(url, params={"data": overpass_query})
                response.raise_for_status()
                data = response.json()

                node_coordinates = {
                    node["id"]: [node["lon"], node["lat"]]
                    for node in data.get("elements", [])
                    if node["type"] == "node"
                }

                
                for element in data.get("elements", []):
                    if feature == "pole":
                        if element["type"] == "node" and "lon" in element and "lat" in element:
                            props = {"name_microgrid": grid_name, "id": element["id"]}
                            tags = element.get("tags")
                            if isinstance(tags, dict):
                                props.update(tags)
                            geom = {"type": "Point","coordinates":[element["lon"],element["lat"]]}
                            geojson_features.append(json.dumps({"type":"Feature","properties":props,"geometry":geom}))
                    else:
                        if element["type"] == "way" and "nodes" in element:
                            coords = [node_coordinates[nid] for nid in element["nodes"] if nid in node_coordinates]
                            if not coords: continue
                            props = {"name_microgrid": grid_name, "id": element["id"]}
                            tags = element.get("tags")
                            if isinstance(tags, dict):
                                props.update(tags)
                            if geometry_type == "Polygon":
                                geom = {"type":"Polygon","coordinates":[coords]}
                            else:
                                geom = {"type":"LineString","coordinates":coords}
                            geojson_features.append(json.dumps({"type":"Feature","properties":props,"geometry":geom}))

            except Exception as e:
                logger.error(f"Request/JSON error for microgrid: {grid_name} feature {feature}: {e}")

        # write output file for this feature
        outpath = Path(path) / filename
        outpath.parent.mkdir(parents=True, exist_ok=True)
        with open(outpath, "w") as f:
            f.write('{"type":"FeatureCollection","features":[\n')
            f.write(",\n".join(geojson_features))
            f.write("\n]}\n")
        logger.info(f"GeoJSON saved to {outpath}")


def download_and_merge_Microsoft_buildings(url, microgrid_list, osm_path, output_path):
    """
    Download and merge Microsoft Building Footprints with existing OSM buildings
    for one or more microgrid bounding boxes.

    The function:
    • reads the Microsoft Buildings tile index CSV from `url` and, for each
        microgrid’s bounding box, determines the quadkeys (zoom level 9) covering it;
    • downloads the corresponding Microsoft building footprints (GeoJSON lines),
        converts them to geometries, filters them to the microgrid extent, and
        concatenates them into one GeoDataFrame with IDs and height attributes;
    • loads an existing OSM buildings file from `osm_path`;
    • reprojects both datasets to a metric CRS, buffers the Microsoft buildings
        slightly, and performs a spatial join with OSM to detect overlaps;
    • cleans and merges attributes from overlapping features with the remaining
        unique OSM buildings;
    • outputs the merged building set as a single GeoJSON FeatureCollection
        written to `output_path`.

    This produces one unified building dataset combining Microsoft and OSM data
    within the given microgrid areas.
    """

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
    cluster_unique = osm_gdf.drop(index=duplicate_idxs).to_crs(crs_meters)

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
    with open(output_path, "w", encoding="utf-8") as f:
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

        for f in out_formats:
            # SALVA AL PLURALE
            new_file_name = Path.joinpath(store_path_resources, f"all_raw_buildings.{f}")
            # accetta sia *buildings.* sia *building.* prodotti da earth_osm
            old_file = (list(Path(out_path).glob(f"*buildings.{f}"))
                        or list(Path(out_path).glob(f"*building.{f}")))

            if not old_file:
                # crea file vuoto valido
                with open(new_file_name, "w", encoding="utf-8") as fp:
                    if f == "geojson":
                        fp.write('{"type":"FeatureCollection","features":[]}\n')
                    else:
                        fp.write("")
            else:
                logger.info(f"Move {old_file[0]} to {new_file_name}")
                shutil.move(old_file[0], new_file_name)

    elif snakemake.config["enable"]["download_osm_method"] == "overpass":
        microgrids_list = snakemake.config["microgrids_list"]
        features = [
            "building",
            "minor_line",
            "cable",
            "generator",
            "substation",
            "pole",
        ]
        overpass_url = "https://overpass-api.de/api/interpreter"
        output_dir = Path.cwd() / "resources" / RDIR / "osm" / "raw"
        retrieve_osm_data_geojson(microgrids_list, features, overpass_url, output_dir)

        if snakemake.config["enable"]["download_and_merge_microsoft_ML_building"]:
            # USA SEMPRE il PLURALE
            osm_path = output_dir / "all_raw_buildings.geojson"
            export_path = output_dir / "all_raw_buildings.geojson"
            microsoft_data_url = "https://minedbuildings.z5.web.core.windows.net/global-buildings/dataset-links.csv"
            download_and_merge_Microsoft_buildings(
                microsoft_data_url, microgrids_list, osm_path, export_path
            )

# -*- coding: utf-8 -*-
import json
import logging
import os
import shutil
import time
from pathlib import Path

import geopandas as gpd
import mercantile
import pandas as pd
import requests
import yaml
from _helpers_dist import configure_logging, create_logger, read_osm_config
from earth_osm import eo
from shapely import geometry
from shapely.geometry import LineString, Point, Polygon, mapping

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


def retrieve_osm_data_geojson(
    microgrids_list,
    url,
    output_path_buildings,
    output_path_minor_lines,
    output_path_cables,
    output_path_generators,
    output_path_substations,
    output_path_poles,
):
    """
    Download OpenStreetMap data via Overpass for a set of microgrid bounding boxes
    and export each requested feature as GeoJSON. For each feature within a given
    bounding box, an Overpass query is built and executed. Geometries are rebuilt
    from returned nodes, and tag columns are reduced to the subset of interest.
    If no data are collected, an empty file is saved and a warning is logged.

    Parameters
    ----------
    microgrids_list : dict
        Dictionary with the bounding-box coordinates for each microgrid.
    url : str
        Overpass API endpoint.
    output_path_buildings : str
        Output path for buildings.
    output_path_minor_lines : str
        Output path for minor_line.
    output_path_cables : str
        Output path for cable.
    output_path_generators : str
        Output path for generator.
    output_path_substations : str
        Output path for substation.
    output_path_poles : str
        Output path for pole.

    Notes
    -----
    Overpass is rate-limited; consider increasing the fixed delay or using mirror
    endpoints if you frequently encounter 429/504 responses.
    """

    def _run_feature(feature, geometry_type, outpath):
        DELAY_S = 2.0
        EMPTY_RETRY_SLEEP_S = 5.0
        last_err: Exception | None = None

        def _collect_records() -> list[dict]:
            records: list[dict] = []
            for grid_name, grid_data in microgrids_list.items():
                lat_min = grid_data["lat_min"]
                lon_min = grid_data["lon_min"]
                lat_max = grid_data["lat_max"]
                lon_max = grid_data["lon_max"]

                if feature == "building":
                    overpass_query = f"""
                    [out:json][timeout:300];
                    ( way["building"]({lat_min},{lon_min},{lat_max},{lon_max}); );
                    (._;>;);
                    out body;"""
                elif feature == "minor_line":
                    overpass_query = f"""
                    [out:json][timeout:300];
                    ( way["power"~"^(minor_line)$"]({lat_min},{lon_min},{lat_max},{lon_max});
                      relation["power"~"^(minor_line)$"]({lat_min},{lon_min},{lat_max},{lon_max}); );
                    (._;>;);
                    out body;"""
                elif feature == "cable":
                    overpass_query = f"""
                    [out:json][timeout:300];
                    ( way["power"~"^(cable)$"]({lat_min},{lon_min},{lat_max},{lon_max});
                      relation["power"~"^(cable)$"]({lat_min},{lon_min},{lat_max},{lon_max}); );
                    (._;>;);
                    out body;"""
                elif feature == "generator":
                    overpass_query = f"""
                    [out:json][timeout:300];
                    ( node["power"~"generator|plant"]({lat_min},{lon_min},{lat_max},{lon_max});
                      way["power"~"generator|plant"]({lat_min},{lon_min},{lat_max},{lon_max});
                      relation["power"~"generator|plant"]({lat_min},{lon_min},{lat_max},{lon_max}); );
                    (._;>;);
                    out body;"""
                elif feature == "substation":
                    overpass_query = f"""
                    [out:json][timeout:300];
                    ( node["power"="substation"]({lat_min},{lon_min},{lat_max},{lon_max});
                      way["power"="substation"]({lat_min},{lon_min},{lat_max},{lon_max}); );
                    (._;>;);
                    out body;"""
                elif feature == "pole":
                    overpass_query = f"""
                    [out:json][timeout:300];
                    ( node["power"="pole"]({lat_min},{lon_min},{lat_max},{lon_max}); );
                    out body;"""
                else:
                    return records

                try:
                    logger.info(f"Overpass query: {grid_name} / {feature}")
                    r = requests.get(
                        url, params={"data": overpass_query}, timeout=240
                    )  # ↑ timeout locale
                    if r.status_code in (429, 502, 503, 504):  # retry minimale
                        time.sleep(5)
                        r = requests.get(
                            url, params={"data": overpass_query}, timeout=240
                        )
                    r.raise_for_status()
                    data = r.json()

                    node_coords = {
                        n["id"]: (n["lon"], n["lat"])
                        for n in data.get("elements", [])
                        if n.get("type") == "node" and "lon" in n and "lat" in n
                    }

                    for el in data.get("elements", []):
                        if (
                            feature == "pole"
                            and el.get("type") == "node"
                            and "lon" in el
                            and "lat" in el
                        ):
                            props = {"name_microgrid": grid_name, "id": el.get("id")}
                            props.update(el.get("tags") or {})
                            records.append(
                                {**props, "geometry": Point(el["lon"], el["lat"])}
                            )
                            continue

                        if el.get("type") == "way" and "nodes" in el:
                            coords = [
                                node_coords[nid]
                                for nid in el["nodes"]
                                if nid in node_coords
                            ]
                            if geometry_type == "Polygon":
                                if len(coords) >= 3:
                                    if coords[0] != coords[-1]:
                                        coords = coords + [coords[0]]
                                    geom = Polygon(coords)
                                    props = {
                                        "name_microgrid": grid_name,
                                        "id": el.get("id"),
                                    }
                                    props.update(el.get("tags") or {})
                                    records.append({**props, "geometry": geom})
                            else:
                                if len(coords) >= 2:
                                    geom = LineString(coords)
                                    props = {
                                        "name_microgrid": grid_name,
                                        "id": el.get("id"),
                                    }
                                    props.update(el.get("tags") or {})
                                    records.append({**props, "geometry": geom})

                except Exception as e:
                    nonlocal last_err
                    last_err = e
                    logger.warning(f"Overpass error for {grid_name}/{feature}: {e}.")
                finally:
                    time.sleep(DELAY_S)
            return records

        records = _collect_records()
        if not records:
            logger.warning(
                f"No data for '{feature}'. Retrying once in {EMPTY_RETRY_SLEEP_S:.0f}s..."
            )
            time.sleep(EMPTY_RETRY_SLEEP_S)
            records = _collect_records()

        outpath = Path(outpath)
        outpath.parent.mkdir(parents=True, exist_ok=True)

        if not records:
            with outpath.open("w", encoding="utf-8") as f:
                json.dump({"type": "FeatureCollection", "features": []}, f)
            msg = f"Empty GeoJSON for '{feature}'. Cause: {last_err or 'no data returned'}"
            logger.warning(msg)
            return

        gdf = gpd.GeoDataFrame(records, geometry="geometry", crs="EPSG:4326")
        keep = [
            "name_microgrid",
            "building",
            "type",
            "height_right",
            "id",
            "voltage",
            "power",
            "circuits",
            "wires",
            "cables",
            "frequency",
            "location",
            "operator",
        ]
        cols = [c for c in keep if c in gdf.columns] + ["geometry"]
        gdf = gdf.loc[:, cols].dropna(axis=1, how="all")
        gdf.to_file(outpath, driver="GeoJSON")
        logger.info(f"GeoJSON saved to {outpath} ({len(gdf)} features)")

    if output_path_buildings:
        _run_feature("building", "Polygon", output_path_buildings)
    if output_path_minor_lines:
        _run_feature("minor_line", "LineString", output_path_minor_lines)
    if output_path_cables:
        _run_feature("cable", "LineString", output_path_cables)
    if output_path_generators:
        _run_feature("generator", "Polygon", output_path_generators)
    if output_path_substations:
        _run_feature("substation", "Polygon", output_path_substations)
    if output_path_poles:
        _run_feature("pole", "Point", output_path_poles)


def download_and_merge_Microsoft_buildings(url, microgrid_list, osm_path, output_path):
    """
    Download and merge Microsoft Building data with existing OSM buildings
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

    Parameters
    ----------
    url : str
        HTTP(S) URL to the Microsoft Buildings tile index CSV
    microgrids_list : dict
        Dictionary with the bounding-box coordinates for each microgrid.
    osm_path : str
        Path to the existing OSM buildings GeoJSON.
    output_path : str
        Output path for the merged buildings GeoJSON.

    Note
    ----------------------
    - Uses zoom level 9 quadkeys to fetch Microsoft tiles covering each box.
    - Projects to a metric CRS for buffering and spatial join,
    then outputs in EPSG:4326.
    - Preserves a minimal set of attributes when available: ["name_microgrid",
    "building", "type", "height_right"].
    """
    link = pd.read_csv(url, dtype=str)
    mML_gdf = gpd.GeoDataFrame()
    idx = 0

    # Iterate microgrid bounding boxes
    for gridname, grid_data in microgrid_list.items():
        lon_min, lat_min, lon_max, lat_max = (
            grid_data["lon_min"],
            grid_data["lat_min"],
            grid_data["lon_max"],
            grid_data["lat_max"],
        )
        microgrid_shape = geometry.box(lon_min, lat_min, lon_max, lat_max)

        # Quadkeys at zoom 9 covering the bbox
        quad_keys = list(
            {
                mercantile.quadkey(tile)
                for tile in mercantile.tiles(
                    lon_min, lat_min, lon_max, lat_max, zooms=9
                )
            }
        )

        # Download and filter buildings per quadkey
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

                # Convert to shapely, clip to microgrid, accumulate
                df_json["geometry"] = df_json["geometry"].apply(geometry.shape)
                gdf = gpd.GeoDataFrame(df_json, geometry="geometry", crs="EPSG:4326")
                gdf = gdf[gdf.geometry.within(microgrid_shape)]
                gdf["id"] = range(idx, idx + len(gdf))
                idx += len(gdf)
                mML_gdf = pd.concat([mML_gdf, gdf], ignore_index=True)

    # Exit early if nothing found
    if mML_gdf.empty:
        print("vuoto")
        return

    # Extract height from nested 'properties' dict (fallback to 0)
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

    # Project both to metric CRS
    crs_meters = "EPSG:32632"
    mML_gdf = mML_gdf.to_crs(crs_meters)
    osm_gdf = osm_gdf.to_crs(crs_meters)

    # Light buffer to improve matching
    mML_buffered = mML_gdf.copy()
    mML_buffered["geometry"] = mML_buffered.geometry.buffer(1.0)

    # Spatial join (right columns get _right suffix, e.g. height_right)
    joined = gpd.sjoin(osm_gdf, mML_buffered, how="inner", predicate="intersects")
    joined = joined.astype({"index_right": int})

    if "height_right" not in joined.columns and "height" in joined.columns:
        joined = joined.rename(columns={"height": "height_right"})

    # Clean attributes from matches (drop geometry/properties duplicates)
    clean_join_rows = []
    for _, row in joined.iterrows():
        geom_ = row["geometry"]
        cleaned_props = {
            k: v
            for k, v in row.items()
            if pd.notna(v) and v is not None and k not in ["geometry", "properties"]
        }
        clean_join_rows.append({**cleaned_props, "geometry": geom_})

    cleaned_join_gdf = gpd.GeoDataFrame(clean_join_rows, crs=crs_meters)

    # Remove duplicates from original OSM and merge back
    duplicate_idxs = joined.index
    cluster_unique = osm_gdf.drop(index=duplicate_idxs).to_crs(crs_meters)
    final_gdf = pd.concat([cleaned_join_gdf, cluster_unique], ignore_index=True)
    # Reproject to WGS84
    final_gdf = final_gdf.to_crs("EPSG:4326")

    if "height_right" not in final_gdf.columns and "height" in final_gdf.columns:
        final_gdf = final_gdf.rename(columns={"height": "height_right"})

    # Keep only a minimal set of columns if present
    if not final_gdf.empty:
        keep = ["name_microgrid", "building", "type", "height_right"]
        cols = [c for c in keep if c in final_gdf.columns]
        cols.append("geometry")
        final_gdf = final_gdf.loc[:, cols].dropna(axis=1, how="all")

    # Save with GeoPandas
    outpath = Path(output_path)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    final_gdf.to_file(outpath, driver="GeoJSON")


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
            new_file_name = Path.joinpath(
                store_path_resources, f"all_raw_buildings.{f}"
            )
            old_file = list(Path(out_path).glob(f"*buildings.{f}")) or list(
                Path(out_path).glob(f"*building.{f}")
            )

            if not old_file:
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

        retrieve_osm_data_geojson(
            microgrids_list,
            overpass_url,
            output_path_buildings=snakemake.output["buildings_resources"],
            output_path_minor_lines=snakemake.output["lines_resources"],
            output_path_cables=snakemake.output["cables_resources"],
            output_path_generators=snakemake.output["generators_resources"],
            output_path_substations=snakemake.output["substations_resources"],
            output_path_poles=snakemake.output["poles_resources"],
        )

        if snakemake.config["enable"]["download_and_merge_microsoft_ML_building"]:
            osm_path = output_dir / "all_raw_buildings.geojson"
            export_path = output_dir / "all_raw_buildings.geojson"
            microsoft_data_url = "https://minedbuildings.z5.web.core.windows.net/global-buildings/dataset-links.csv"
            download_and_merge_Microsoft_buildings(
                microsoft_data_url, microgrids_list, osm_path, export_path
            )

# -*- coding: utf-8 -*-
"""
Filter and prepare network data for PyPSA-Distribution.

This script:
1. Renames all buses in the base network using a microgrid prefix.
2. Renames the 'bus' coordinate in all renewable profile .nc files.
3. Detects external grid-connection buses and marks them in the network.
4. Saves the updated network as base_update.nc.

Inputs come from the Snakemake rule `filter_data`.
"""

import logging
import os

import geopandas as gpd
import pypsa
import xarray as xr
from _helpers_dist import configure_logging, sets_path_to_root
from shapely.geometry import Point

logger = logging.getLogger(__name__)


def rename_microgrid_buses(n, microgrids_list):
    """
    Rename all buses in a PyPSA network using a microgrid prefix.

    Parameters
    ----------
    n : pypsa.Network
        The input network.

    microgrids_list : dict
        Dictionary specifying microgrid names. The first key is used
        as the bus prefix.

    Returns
    -------
    n : pypsa.Network
        Network with renamed buses and updated references.
    """
    microgrid_name = list(microgrids_list.keys())[0]

    # Create mapping old_bus → microgrid_1_bus_old_bus
    rename_map = {old: f"{microgrid_name}_bus_{old}" for old in n.buses.index}

    # Rename bus index itself
    n.buses.index = n.buses.index.to_series().map(rename_map)

    # Update all components that contain bus references
    for comp in [
        "lines",
        "links",
        "transformers",
        "loads",
        "generators",
        "stores",
        "storage_units",
    ]:
        if hasattr(n, comp):
            df = getattr(n, comp)
            for col in df.columns:
                if col.startswith("bus"):
                    df[col] = df[col].replace(rename_map)

    logger.info(f"Renamed all buses using prefix '{microgrid_name}_bus_'.")
    return n


def rename_profile_buses(nc_path, microgrid_name):
    """
    Rename the 'bus' dimension of a renewable profile NetCDF file and
    safely overwrite the original file using a temporary write + os.replace.
    """

    ds = xr.open_dataset(nc_path)

    # Skip if no bus dimension
    if "bus" not in ds.coords:
        logger.warning(f"File {nc_path} has no 'bus' dimension — skipping.")
        return

    old_buses = ds.bus.values.astype(str)
    rename_map = {b: f"{microgrid_name}_bus_{b}" for b in old_buses}

    # Apply renaming
    ds = ds.assign_coords(bus=[rename_map[b] for b in old_buses])

    # use file temp to avoid problem of permission denied
    tmp_path = nc_path + ".tmp"

    # Write temporary NC file inside WSL
    ds.to_netcdf(tmp_path)

    # Replace original file (Windows accepts atomic replace)
    import os

    os.replace(tmp_path, nc_path)

    logger.info(f"Updated bus coordinates in {nc_path}.")
    return rename_map


def find_external_connection_buses(n, lines_path, shape_path, crs="EPSG:4326"):
    """
    Identify buses inside the microgrid that connect to lines crossing
    outside the specified shape area.

    Parameters
    ----------
    n : pypsa.Network
        Network containing buses with coordinates.

    lines_path : str
        Path to raw OSM line geometries.

    shape_path : str
        Path to the microgrid boundary polygon.

    crs : str
        CRS for geometry comparison.

    Returns
    -------
    set
        Set of bus indices connecting to external lines.
    """
    lines_raw = gpd.read_file(lines_path).to_crs(crs)
    shape = gpd.read_file(shape_path).to_crs(crs)
    area = shape.unary_union

    # Identify lines that leave the microgrid boundary
    inside_mask = lines_raw.geometry.within(area)
    outside_lines = lines_raw.loc[~inside_mask].copy()

    # Compute line endpoints
    outside_lines["start"] = outside_lines.geometry.apply(lambda g: Point(g.coords[0]))
    outside_lines["end"] = outside_lines.geometry.apply(lambda g: Point(g.coords[-1]))
    outside_lines["start_inside"] = outside_lines["start"].apply(area.contains)
    outside_lines["end_inside"] = outside_lines["end"].apply(area.contains)

    # Convert buses to GeoDataFrame
    buses = gpd.GeoDataFrame(
        n.buses.copy(),
        geometry=gpd.points_from_xy(n.buses.x, n.buses.y),
        crs=crs,
    )

    # Helper: nearest bus to a point
    def nearest_bus(point):
        buses_utm = buses.to_crs("EPSG:3857")
        p_utm = gpd.GeoSeries([point], crs=crs).to_crs("EPSG:3857").iloc[0]
        return buses_utm.geometry.distance(p_utm).idxmin()

    # Collect inside buses that touch outside lines
    external_buses = set()
    for _, row in outside_lines.iterrows():
        if row["start_inside"]:
            external_buses.add(nearest_bus(row["start"]))
        elif row["end_inside"]:
            external_buses.add(nearest_bus(row["end"]))

    return external_buses


def mark_external_buses(n, lines_path, shape_path):
    """
    Mark buses that connect to lines outside the microgrid boundary by
    adding a boolean column n.buses["outside"].

    Parameters
    ----------
    n : pypsa.Network

    lines_path : str
        OSM cleaned lines.

    shape_path : str
        Microgrid area polygon.

    Returns
    -------
    n : pypsa.Network
        Updated network.
    """
    external_ids = find_external_connection_buses(n, lines_path, shape_path)

    n.buses["outside"] = n.buses.index.isin(external_ids)
    logger.info(f"Identified {len(external_ids)} external connection buses.")

    return n


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_dist import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("filter_data")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)

    # Load base network
    n_base = pypsa.Network(snakemake.input.base_network)

    # Extract microgrid name
    microgrid_name = list(snakemake.config["microgrids_list"].keys())[0]

    # 1) Rename buses in the PyPSA network
    n = rename_microgrid_buses(n_base, snakemake.config["microgrids_list"])

    # 2) Rename buses inside all renewable profile .nc files
    for key, path in snakemake.input.items():
        if key.startswith("profile_"):
            rename_profile_buses(path, microgrid_name)

    # 3) Mark external buses
    n = mark_external_buses(n, snakemake.input["raw_lines"], snakemake.input["shape"])

    # 4) Save final updated network
    n.export_to_netcdf(snakemake.output["base_update"])
    logger.info(f"Saved updated network to {snakemake.output['base_update']}")

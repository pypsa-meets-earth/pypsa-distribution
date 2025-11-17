# -*- coding: utf-8 -*-
import logging
import os

import geopandas as gpd
import pypsa
from _helpers_dist import configure_logging, sets_path_to_root
from shapely.geometry import Point

logger = logging.getLogger(__name__)


def rename_microgrid_buses(n, microgrids_list, output_path):
    """
    Rename all buses in the PyPSA network with the given microgrid prefix.
    Parameters
    ----------
    n : pypsa.Network
        Input network.
    microgrid_name : str
        Prefix to add to each bus (e.g., 'microgrid_1').
    Returns
    -------
    n : pypsa.Network
        Updated network with renamed buses and coherent references.
    """
    microgrid_name = list(microgrids_list.keys())[0]
    rename_map = {old: f"{microgrid_name}_bus_{old}" for old in n.buses.index}
    # Rename bus index
    n.buses.index = n.buses.index.to_series().map(rename_map)
    # Update bus references in all components
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
            setattr(n, comp, df)
    # Check integrity
    if hasattr(n, "lines"):
        assert all(n.lines.bus0.isin(n.buses.index))
        assert all(n.lines.bus1.isin(n.buses.index))
    logger.info(f"Buses renamed with prefix '{microgrid_name}_bus_'.")
    # Export updated network
    n.export_to_netcdf(output_path)
    logger.info(f"Network exported to {output_path}")
    return n


def find_external_connection_buses(n, lines_path, shape_path, crs="EPSG:4326"):
    """
    Identify internal buses connected to lines that exit a given area (shape file).
    """
    # Load files and unify CRS
    lines_raw = gpd.read_file(lines_path).to_crs(crs)
    shape = gpd.read_file(shape_path).to_crs(crs)
    area = shape.unary_union
    # Find lines that cross or exit the shape area
    inside_mask = lines_raw.geometry.within(area)
    outside_raw_lines = lines_raw.loc[~inside_mask].copy()
    # Compute line endpoints
    outside_raw_lines["start"] = outside_raw_lines.geometry.apply(
        lambda g: Point(g.coords[0])
    )
    outside_raw_lines["end"] = outside_raw_lines.geometry.apply(
        lambda g: Point(g.coords[-1])
    )
    # Identify which endpoint is inside the area
    outside_raw_lines["start_inside"] = outside_raw_lines["start"].apply(
        lambda p: area.contains(p)
    )
    outside_raw_lines["end_inside"] = outside_raw_lines["end"].apply(
        lambda p: area.contains(p)
    )
    # Prepare buses GeoDataFrame
    buses = n.buses.copy()
    buses = gpd.GeoDataFrame(
        buses, geometry=gpd.points_from_xy(buses.x, buses.y), crs=crs
    )

    # Helper function to find the nearest bus to a given point
    def nearest_bus(point, buses):
        buses_utm = buses.to_crs("EPSG:3857")
        p_utm = gpd.GeoSeries([point], crs=crs).to_crs("EPSG:3857").iloc[0]
        distances = buses_utm.geometry.distance(p_utm)
        return distances.idxmin()

    # Find internal buses connected to external lines
    bus_ids_outside = set()
    for _, row in outside_raw_lines.iterrows():
        if row["start_inside"]:
            p = row["start"]
        elif row["end_inside"]:
            p = row["end"]
        else:
            continue
        bus_ids_outside.add(nearest_bus(p, buses))
    return bus_ids_outside


def mark_external_buses(n, lines_path, shape_path, output_path, crs="EPSG:4326"):
    """
    Identify external connection buses and add a boolean attribute 'outside'
    to the buses of a PyPSA network. The updated network is then saved to disk.

    This function internally calls `find_external_connection_buses` to detect
    which buses are connected to lines that exit the defined area.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network to modify.
    lines_path : str
        Path to the GeoJSON/Shapefile containing line geometries.
    shape_path : str
        Path to the GeoJSON/Shapefile defining the network area.
    output_path : str
        Path to export the updated network (.nc file).
    crs : str, optional
        Coordinate Reference System (default: "EPSG:4326").

    Returns
    -------
    n : pypsa.Network
        Updated network with an 'outside' column added to n.buses.
    """
    # Identify buses connected to external lines
    bus_ids_outside = find_external_connection_buses(n, lines_path, shape_path, crs)
    # Add 'outside' flag to buses
    n.buses["outside"] = n.buses.index.isin(bus_ids_outside)
    logger.info(
        f"Added 'outside' column to n.buses: {n.buses['outside'].sum()} external buses identified."
    )
    # Export updated network
    n.export_to_netcdf(output_path)
    logger.info(f"Updated network saved to {output_path}")
    return n


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_dist import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("filter_data")
        sets_path_to_root("pypsa-distribution")

    # Configura logging standard del workflow
    configure_logging(snakemake)

    n_base = pypsa.Network(snakemake.input.base_network)
    output = snakemake.output["base_update"]
    n = rename_microgrid_buses(
        n_base,
        snakemake.config["microgrids_list"],
        output,
    )

    mark_external_buses(
        n, snakemake.input["raw_lines"], snakemake.input["shape"], output
    )

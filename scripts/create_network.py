# -*- coding: utf-8 -*-

import json
import logging
import os
from itertools import combinations

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa
from _helpers_dist import configure_logging, read_geojson, sets_path_to_root
from geopy import distance as geodistance
from scipy.spatial import Delaunay, distance
from shapely.geometry import Polygon

_logger = logging.getLogger(__name__)
_logger.setLevel(logging.INFO)


def create_network():
    """
    Creates a PyPSA network and sets the snapshots for the network
    """
    n = pypsa.Network()

    # Set the name of the network
    n.name = "PyPSA-Distribution"

    # Set the snapshots for the network
    n.set_snapshots(pd.date_range(freq="h", **snakemake.config["snapshots"]))

    # Normalize the snapshot weightings
    n.snapshot_weightings[:] *= 8760.0 / n.snapshot_weightings.sum()

    # Return the created network
    return n


def create_microgrid_network(
    n, input_file, number_microgrids, voltage_level, line_type
):
    """
    Creates local microgrid networks within the PyPSA network. The local microgrid networks are distribution networks created based on
    the buildings data, stored in "resources/buildings/microgrids_buildings.geojson". Then the buses are connected together through lines
    according to the output of a Delaunay Triangulation.
    """
    # Load the GeoJSON file

    with open(input_file) as f:
        data = json.load(f)

    # Keep track of the bus coordinates and microgrid IDs
    bus_coords = set()
    number_microgrids = len(number_microgrids.keys())
    microgrid_ids = [f"microgrid_{i+1}" for i in range(number_microgrids)]
    # microgrid_ids = set()

    # Iterate over each feature in the GeoDataFrame
    for feature in data["features"]:
        # Get the point geometry
        point_geom = feature["geometry"]

        # Create a bus at the point location with microgrid ID included in bus name
        bus_name = f"bus_{feature['properties']['cluster']}"

        x, y = point_geom["coordinates"][0], point_geom["coordinates"][1]

        # Check for overlapping microgrids and raise an error if happening
        if (x, y) in bus_coords:
            raise ValueError(
                "Overlapping microgrids detected, adjust the coordinates in the config.yaml file"
            )

        # Add the buses to the network and update the set of bus coordinates and microgrid IDs
        n.add("Bus", bus_name, x=x, y=y, v_nom=voltage_level)
        bus_coords.add((x, y))

    # Iterate over each microgrid
    for microgrid_id in microgrid_ids:
        coords = np.column_stack((n.buses.x.values, n.buses.y.values))

        # Create a Delaunay triangulation of the bus coordinates
        tri = Delaunay(coords)

        # Remove edges that connect the same pair of buses
        edges = []
        for simplex in tri.simplices:
            for i in range(3):
                if i < 2:
                    edge = sorted([simplex[i], simplex[i + 1]])
            else:
                edge = sorted([simplex[i], simplex[0]])
            if edge not in edges:
                edges.append(edge)

        # # Create a matrix of bus coordinates

        # # Create a Delaunay triangulation of the bus coordinates
        # tri = Delaunay(coords)
        # edges = tri.simplices[(tri.simplices < len(coords)).all(axis=1)]

        line_type = line_type

    line_name = [f"{microgrid_id}_line_{i}_{j}" for (i, j) in edges]

    edges_list = [(i, j) for (i, j) in edges]
    i_bus0 = [k[0] for k in edges_list]
    i_bus1 = [k[1] for k in edges_list]
    # keep order of buses to evaluate length of the lines
    lines_df = pd.concat(
        [
            n.buses.loc[n.buses.index[i_bus0], ["x", "y"]].reset_index(),
            n.buses.loc[n.buses.index[i_bus1], ["x", "y"]].reset_index(),
        ],
        axis=1,
        ignore_index=True,
    ).set_axis(["bus0", "x0", "y0", "bus1", "x1", "y1"], axis=1)

    lines_df["length"] = lines_df.apply(
        lambda x: geodistance.geodesic((x["x0"], x["y0"]), (x["x1"], x["y1"])).meters
        / 1e3,
        axis=1,
    )

    n.madd(
        "Line",
        line_name,
        bus0=lines_df["bus0"].values,
        bus1=lines_df["bus1"].values,
        type=line_type,
        length=lines_df["length"],
        # Default settings imply zero s_nom and non-extandable lines
        s_nom=0.1,  # quick fix
        # s_nom_extendable=True # makes optimisation longer
    )


def add_bus_at_center(n, number_microgrids, voltage_level, line_type):
    """
    Adds a new bus to each network at the center of the existing buses.
    This is the bus to which the generation, the storage and the load will be attached.
    """
    number_microgrids = len(number_microgrids.keys())
    microgrid_ids = [f"microgrid_{i+1}" for i in range(number_microgrids)]

    # Iterate over each microgrid
    for microgrid_id in microgrid_ids:
        # Select the buses belonging to this microgrid
        microgrid_buses = n.buses.loc[
            n.buses.index.str.startswith(f"{microgrid_id}_bus_")
        ]

        # Create a matrix of bus coordinates
        coords = np.column_stack((microgrid_buses.x.values, microgrid_buses.y.values))
        polygon = Polygon(coords)
        s = gpd.GeoSeries(polygon)
        s = s.centroid

        # Create a new bus at the centroid
        center_bus_name = f"new_bus_{microgrid_id}"
        n.add(
            "Bus",
            center_bus_name,
            x=float(s.x.iloc[0]),
            y=float(s.y.iloc[0]),
            v_nom=voltage_level,
        )

        # Find the two closest buses to the new bus
        closest_buses = microgrid_buses.iloc[
            distance.cdist([(float(s.x.iloc[0]), float(s.y.iloc[0]))], coords).argmin()
        ]
        closest_buses = closest_buses.iloc[[0, 1]]
        line_type = line_type

        # Add lines to connect the new bus to the closest buses)

        # Add lines to connect the new bus to the closest buses
        for _, bus in closest_buses.to_frame().iterrows():
            line_name = f"{microgrid_id}_line_{center_bus_name}_{bus.name}"
            x1, y1 = n.buses.loc[bus.index].x, n.buses.loc[bus.index].y
            x2, y2 = n.buses.loc[center_bus_name].x, n.buses.loc[center_bus_name].y
            length = ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5
            n.add(
                "Line",
                line_name,
                bus0=center_bus_name,
                bus1=bus.index,
                type=line_type,
                length=length,
            )


def plot_microgrid_network(n):
    # Create a new figure and axis
    fig, ax = plt.subplots()

    # Plot each bus in the network
    for bus_name, bus in n.buses.iterrows():
        ax.plot(bus.x, bus.y, "o", color="blue")

    # Plot each line in the network
    for line_name, line in n.lines.iterrows():
        bus0 = n.buses.loc[line.bus0]
        bus1 = n.buses.loc[line.bus1]
        ax.plot([bus0.x, bus1.x], [bus0.y, bus1.y], "-", color="black")

    # Set the axis limits to include all buses in the network
    ax.set_xlim(n.buses.x.min() - 0.1, n.buses.x.max() + 0.1)
    ax.set_ylim(n.buses.y.min() - 0.1, n.buses.y.max() + 0.1)

    # Set the title and labels for the plot
    ax.set_title("Networks of the microgrids")
    ax.set_xlabel("X Coordinate")
    ax.set_ylabel("Y Coordinate")

    # Show the plot
    plt.show()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_dist import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("create_network")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)

    n = create_network()

    create_microgrid_network(
        n,
        snakemake.input["clusters"],
        snakemake.config["microgrids_list"],
        snakemake.config["electricity"]["voltage"],
        snakemake.config["electricity"]["line_type"],
    )

    # add_bus_at_center(n,
    #                   snakemake.config["microgrids_list"],
    #                   snakemake.config["electricity"]["voltage"],
    #                   snakemake.config["electricity"]["line_type"])

    # plot_microgrid_network(n)
    a = 12
    n.export_to_netcdf(snakemake.output[0])

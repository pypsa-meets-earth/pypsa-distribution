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
from pyproj import Transformer
from scipy.spatial import Delaunay, distance
from shapely.geometry import Point, Polygon

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


def calculate_power_node_position(load_file, cluster_bus):
    load_sums = load_file.sum(numeric_only=True)
    load_sums.index = cluster_bus["cluster"]
    gdf = cluster_bus.set_index("cluster")
    gdf["cluster_load"] = load_sums.values.T
    x_wgt_avg = (gdf.geometry.x * load_sums).sum() / load_sums.sum()
    y_wgt_avg = (gdf.geometry.y * load_sums).sum() / load_sums.sum()

    return x_wgt_avg, y_wgt_avg


def create_microgrid_network(
    n, input_file, voltage_level, line_type, microgrid_list, input_path
):
    """
    Creates local microgrid networks within the PyPSA network. The local microgrid networks are distribution networks created based on
    the buildings data, stored in "resources/buildings/microgrids_buildings.geojson".
    Each bus corresponds to a cluster of buildings within a microgrid, with its coordinates defined in the input GeoJSON file.
    The lines connecting buses are determined using Delaunay triangulation, ensuring minimal total line length.
    The function avoids duplicate buses and ensures buses are assigned to the correct SubNetwork.
        Parameters
    ----------
    n : pypsa.Network
        The PyPSA network object to which microgrid buses and lines will be added.
    input_file : str
        Path to the GeoJSON file containing building and microgrid data.
    voltage_level : float
        The nominal voltage level to be assigned to the buses.
    line_type : str
        The type of lines to be used for connecting the buses (e.g., "AC").
    microgrid_list : dict
        A dictionary containing the list of microgrids. Keys are microgrid names,
        and values are metadata about each microgrid.
    Output
    ------
    The PyPSA network (`n`) is updated with:
    - Buses for each microgrid, identified by cluster ID and associated with a SubNetwork.
    - Lines connecting buses within each microgrid based on Delaunay triangulation.
    """

    data = gpd.read_file(input_file)
    load = pd.read_csv(input_path)
    bus_coords = set()  # Keep track of bus coordinates to avoid duplicates

    for grid_name, grid_data in microgrid_list.items():
        # List to store bus names and their positions for triangulation
        microgrid_buses = []
        bus_positions = []

        # Filter data for the current microgrid
        grid_data = data[data["name_microgrid"] == grid_name]
        load_data = load[[col for col in load.columns if grid_name in col]]
        x_gen_bus, y_gen_bus = calculate_power_node_position(load_data, grid_data)
        gen_bus_name = f"{grid_name}_gen_bus"
        n.add(
            "Bus",
            gen_bus_name,
            x=x_gen_bus,
            y=y_gen_bus,
            v_nom=voltage_level,
            sub_network=grid_name,
        )
        microgrid_buses.append(gen_bus_name)
        bus_positions.append((x_gen_bus, y_gen_bus))

        # Create a SubNetwork for the current microgrid if it does not exist
        if grid_name not in n.sub_networks.index:
            n.add("SubNetwork", grid_name, carrier="electricity")

        for _, feature in grid_data.iterrows():
            point_geom = feature.geometry
            bus_name = f"{grid_name}_bus_{feature['cluster']}"
            x, y = point_geom.x, point_geom.y

            # Skip duplicate buses or overlapping coordinates
            if bus_name in n.buses.index:
                continue
            if (x, y) in bus_coords:
                raise ValueError(
                    f"Overlapping microgrids detected at {x}, {y}. Adjust the configuration."
                )

            # Add the bus to the network and assign it to the SubNetwork
            n.add("Bus", bus_name, x=x, y=y, v_nom=voltage_level, sub_network=grid_name)
            bus_coords.add((x, y))
            microgrid_buses.append(bus_name)
            bus_positions.append((x, y))

        # Check if there are enough points for triangulation
        if len(bus_positions) < 3:
            print(f"Not enough points for triangulation in {grid_name}.")
            continue

        # Perform Delaunay triangulation to determine bus connections
        coords = np.array(bus_positions)
        tri = Delaunay(coords)

        # Collect unique edges from the Delaunay triangulation
        edges = set()
        for simplex in tri.simplices:
            for i in range(3):
                edge = tuple(sorted([simplex[i], simplex[(i + 1) % 3]]))
                edges.add(edge)

        # Add lines to the network based on the triangulation edges
        for i, j in edges:
            bus0 = microgrid_buses[i]
            bus1 = microgrid_buses[j]
            line_name = f"{grid_name}_line_{bus0}_{bus1}"

            # Skip if the line already exists
            if line_name in n.lines.index:
                continue

            # Retrieve the coordinates of the buses
            x1, y1 = n.buses.loc[bus0].x, n.buses.loc[bus0].y
            x2, y2 = n.buses.loc[bus1].x, n.buses.loc[bus1].y

            # Calculate the distance between buses (in km)
            length = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2) / 1000

            # Add the line to the network
            n.add(
                "Line",
                line_name,
                bus0=bus0,
                bus1=bus1,
                type=line_type,
                length=length,
                s_nom=0.1,
                s_nom_extendable=True,
            )

    if len(microgrid_list.keys()) == 2:
        bus_positions = []
        microgrids = []
        for each in microgrid_list.keys():
            x = n.buses.loc[f"{each}_gen_bus"].x
            y = n.buses.loc[f"{each}_gen_bus"].y
            bus_positions.append((x, y))
            gen_bus_name = f"{each}_gen_bus"
            microgrids.append(gen_bus_name)
        bus0 = microgrids[0]
        bus1 = microgrids[1]
        line_name = f"interconnection_line_between_{bus0}_and_{bus1}"
        x1, y1 = n.buses.loc[bus0].x, n.buses.loc[bus0].y
        x2, y2 = n.buses.loc[bus1].x, n.buses.loc[bus1].y
        length = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2) / 1000
        n.add(
                "Line",
                line_name,
                bus0=bus0,
                bus1=bus1,
                type=line_type,
                length=length,
                s_nom=0.1,
                s_nom_extendable=True,
                )
    elif len(microgrid_list.keys()) > 2:
                bus_positions = []
                microgrids = []
                for each in microgrid_list.keys():
                    x = n.buses.loc[f"{each}_gen_bus"].x
                    y = n.buses.loc[f"{each}_gen_bus"].y
                    bus_positions.append((x, y))
                    gen_bus_name = f"{each}_gen_bus"
                    microgrids.append(gen_bus_name)
                # add lines between generation buses
                coords = np.array(bus_positions)
                tri = Delaunay(coords)
                # Collect unique edges from the Delaunay triangulation
                edges = set()
                for simplex in tri.simplices:
                    for i in range(3):
                        edge = tuple(sorted([simplex[i], simplex[(i + 1) % 3]]))
                        edges.add(edge)
            # Add lines to the network based on the triangulation edges
                for i, j in edges:
                    bus0 = microgrids[i]
                    bus1 = microgrids[j]
                    line_name = f"interconnection_line_between_{bus0}_and_{bus1}"
                    if line_name in n.lines.index:
                        continue
                    x1, y1 = n.buses.loc[bus0].x, n.buses.loc[bus0].y
                    x2, y2 = n.buses.loc[bus1].x, n.buses.loc[bus1].y
                    length = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2) / 1000
                    # Add the line to the network
                    n.add(
                        "Line",
                        line_name,
                        bus0=bus0,
                        bus1=bus1,
                        type=line_type,
                        length=length,
                        s_nom=0.1,
                        s_nom_extendable=True,
                        )


# def add_bus_at_center(n, number_microgrids, voltage_level, line_type):
#     """
#     Adds a new bus to each network at the center of the existing buses.
#     This is the bus to which the generation, the storage and the load will be attached.
#     """
#     number_microgrids = len(number_microgrids.keys())
#     microgrid_ids = [f"microgrid_{i+1}" for i in range(number_microgrids)]

#     # Iterate over each microgrid
#     for microgrid_id in microgrid_ids:
#         # Select the buses belonging to this microgrid
#         microgrid_buses = n.buses.loc[
#             n.buses.index.str.startswith(f"{microgrid_id}_bus_")
#         ]

#         # Create a matrix of bus coordinates
#         coords = np.column_stack((microgrid_buses.x.values, microgrid_buses.y.values))
#         polygon = Polygon(coords)
#         s = gpd.GeoSeries(polygon)
#         s = s.centroid

#         # Create a new bus at the centroid
#         center_bus_name = f"new_bus_{microgrid_id}"
#         n.add(
#             "Bus",
#             center_bus_name,
#             x=float(s.x.iloc[0]),
#             y=float(s.y.iloc[0]),
#             v_nom=voltage_level,
#         )

#         # Find the two closest buses to the new bus
#         closest_buses = microgrid_buses.iloc[
#             distance.cdist([(float(s.x.iloc[0]), float(s.y.iloc[0]))], coords).argmin()
#         ]
#         closest_buses = closest_buses.iloc[[0, 1]]
#         line_type = line_type

#         # Add lines to connect the new bus to the closest buses)

#         # Add lines to connect the new bus to the closest buses
#         for _, bus in closest_buses.to_frame().iterrows():
#             line_name = f"{microgrid_id}_line_{center_bus_name}_{bus.name}"
#             x1, y1 = n.buses.loc[bus.index].x, n.buses.loc[bus.index].y
#             x2, y2 = n.buses.loc[center_bus_name].x, n.buses.loc[center_bus_name].y
#             length = ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5
#             n.add(
#                 "Line",
#                 line_name,
#                 bus0=center_bus_name,
#                 bus1=bus.index,
#                 type=line_type,
#                 length=length,
#             )


# def plot_microgrid_network(n):
#     # Create a new figure and axis
#     fig, ax = plt.subplots()

#     # Plot each bus in the network
#     for bus_name, bus in n.buses.iterrows():
#         ax.plot(bus.x, bus.y, "o", color="blue")

#     # Plot each line in the network
#     for line_name, line in n.lines.iterrows():
#         bus0 = n.buses.loc[line.bus0]
#         bus1 = n.buses.loc[line.bus1]
#         ax.plot([bus0.x, bus1.x], [bus0.y, bus1.y], "-", color="black")

#     # Set the axis limits to include all buses in the network
#     ax.set_xlim(n.buses.x.min() - 0.1, n.buses.x.max() + 0.1)
#     ax.set_ylim(n.buses.y.min() - 0.1, n.buses.y.max() + 0.1)

#     # Set the title and labels for the plot
#     ax.set_title("Networks of the microgrids")
#     ax.set_xlabel("X Coordinate")
#     ax.set_ylabel("Y Coordinate")

#     # Show the plot
#     plt.show()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_dist import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("create_network")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)

    n = create_network()
    microgrids_list = snakemake.config["microgrids_list"]

    create_microgrid_network(
        n,
        snakemake.input["clusters"],
        snakemake.config["electricity"]["voltage"],
        snakemake.config["electricity"]["line_type"],
        microgrids_list,
        snakemake.input["load"],
    )
    a = 12
    n.export_to_netcdf(snakemake.output[0])

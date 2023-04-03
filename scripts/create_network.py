# -*- coding: utf-8 -*-
import json
import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa
from _helpers_dist import configure_logging, sets_path_to_root
from scipy.spatial import Delaunay
from shapely.geometry import shape

_logger = logging.getLogger(__name__)
_logger.setLevel(logging.INFO)


def create_network():
    n = pypsa.Network()

    # Set the name of the network
    n.name = "PyPSA-Distribution"

    # Set the snapshots for the network
    n.set_snapshots(pd.date_range(freq="h", **snakemake.config["snapshots"]))

    # Normalize the snapshot weightings
    n.snapshot_weightings[:] *= 8760.0 / n.snapshot_weightings.sum()

    # Return the created network
    return n


def create_microgrid_network(n, input_file):
    # Load the GeoJSON file
    with open(input_file) as f:
        data = json.load(f)

    # Keep track of the bus coordinates and microgrid IDs
    bus_coords = set()
    microgrid_ids = set()

    # Iterate over each feature in the GeoJSON file
    for feature in data["features"]:
        # Get the point geometry
        point_geom = shape(feature["geometry"])

        # Create a bus at the point location with microgrid ID included in bus name
        bus_name = f"{feature['properties']['microgrid_id']}_bus_{feature['id']}"
        x, y = point_geom.x, point_geom.y

        # Check for overlapping microgrids and raise an error if happening
        if (x, y) in bus_coords:
            raise ValueError(
                "Overlapping microgrids detected, adjust the coordinates in the config.yaml file"
            )

        # Add the bus to the network and update the set of bus coordinates and microgrid IDs
        n.add("Bus", bus_name, x=x, y=y, v_nom=0.220)
        bus_coords.add((x, y))
        microgrid_ids.add(feature["properties"]["microgrid_id"])

    # Iterate over each microgrid
    for microgrid_id in microgrid_ids:
        # Select the buses belonging to this microgrid
        microgrid_buses = n.buses.loc[
            n.buses.index.str.startswith(f"{microgrid_id}_bus_")
        ]

        # Create a matrix of bus coordinates
        coords = np.column_stack((microgrid_buses.x.values, microgrid_buses.y.values))

        # Create a Delaunay triangulation
        tri = Delaunay(coords)

        # Add the edges of the triangulation to the network as lines, but only between buses in the same microgrid
        for edge in tri.simplices:
            if (edge[0] < len(microgrid_buses)) and (edge[1] < len(microgrid_buses)):
                n.add(
                    "Line",
                    f"{microgrid_id}_line_{edge[0]}_{edge[1]}",
                    bus0=microgrid_buses.index[edge[0]],
                    bus1=microgrid_buses.index[edge[1]],
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

    create_microgrid_network(n, 
                            snakemake.input["microgrids_buildings"])
    
    #plot_microgrid_network(n)

    n.export_to_netcdf(snakemake.output[0])

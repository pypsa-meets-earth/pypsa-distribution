# -*- coding: utf-8 -*-
"""
Creates a base network with one bus

Relevant Settings
-----------------
.. code:: yaml
    snapshots:

Inputs
------
Outputs
-------
- ``networks/base_{i}.nc``
   
Description
-----------
This script creates a PyPSA network with one AC bus for each of the microgrids contained in config.yaml
"""

import os

import pandas as pd
import pypsa
from _helpers import configure_logging, sets_path_to_root


def create_networks(num_networks):
    networks = []

    for i in range(num_networks):
        n = pypsa.Network()

        # Set the name of the network
        n.name = f"PyPSA-Distribution-{i+1}"

        # Set the snapshots for the network
        n.set_snapshots(pd.date_range(freq="h", **snakemake.config["snapshots"]))

        # Normalize the snapshot weightings
        n.snapshot_weightings[:] *= 8760.0 / n.snapshot_weightings.sum()

        networks.append(n)

    # Return the list of created networks
    return networks


def add_bus_to_networks(networks):
    # Add one AC bus to each network
    for n in networks:
        n.madd("Bus", ["onebus"], carrier="AC", v_nom=0.220)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("create_network")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)

    number_microgrids = len(snakemake.config["microgrids_list"])
    networks = create_networks(number_microgrids)

    add_bus_to_networks(networks)

    for i, n in enumerate(networks):
        # Export each network to a separate file, using the index i in the filename
        output_filename = f"networks/base_{i+1}.nc"
        n.export_to_netcdf(output_filename)

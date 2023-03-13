# -*- coding: utf-8 -*-
"""
Creates a base network with buses

Relevant Settings
-----------------
.. code:: yaml
    snapshots:
    microgrids_list

Inputs
------
Outputs
-------
- ``networks/base.nc``
   
Description
-----------
This script creates a PyPSA network with as much AC buses as microgrids specified in the file config.yaml
"""

import os

import pandas as pd
import pypsa
from _helpers import configure_logging, sets_path_to_root


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


def add_buses_to_network(n, number_microgrids):
    # Add buses to the network based on the number of microgrids
    for i in range(number_microgrids):
        n.madd("Bus", [f"bus{i+1}"], carrier="AC", v_nom=0.220)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("create_network")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)

    number_microgrids = len(snakemake.config["microgrids_list"])

    n = create_network()
    add_buses_to_network(n, number_microgrids)
    print(n)
    n.export_to_netcdf(snakemake.output[0])

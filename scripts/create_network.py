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
- ``networks/base.nc``
   
Description
-----------
This script creates a PyPSA network with one AC bus.
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


def add_bus_to_network(n):
    # Add one AC bus to the network
    n.madd("Bus", ["onebus"], carrier="AC", v_nom=0.220)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("create_network")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)

    n = create_network()
    add_bus_to_network(n)
    n.export_to_netcdf(snakemake.output[0])

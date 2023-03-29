import json
import pypsa
from scipy.spatial import Delaunay
import numpy as np
import matplotlib as plt
from shapely.geometry import shape
import pandas as pd
import os
from _helpers_dist import (
    configure_logging,
    sets_path_to_root,
)
import logging


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


# Iterate over each feature in the GeoJSON file
    for feature in data['features']:
    # Get the point geometry
        point_geom = shape(feature['geometry'])

    # Create a bus at the point location with microgrid ID included in bus name
        bus_name = f"{feature['properties']['microgrid_id']}_bus_{feature['id']}"
        n.add("Bus", bus_name, x=point_geom.x, y=point_geom.y, v_nom=0.220)

# Group the buses by microgrid ID
    bus_groups = n.buses.groupby(lambda bus: bus.split("_")[0])

# Iterate over each microgrid
    for microgrid_id, buses in bus_groups:
    # Select the buses belonging to this microgrid
        microgrid_buses = n.buses.loc[buses.index[buses.index.str.contains(microgrid_id)]]

    # Create a matrix of bus coordinates
        coords = np.column_stack((microgrid_buses.x.values, microgrid_buses.y.values))

    # Create a Delaunay triangulation
        tri = Delaunay(coords)

    # Add the edges of the triangulation to the network as lines
        for edge in tri.simplices:
            n.add("Line", f"{microgrid_id}_line_{edge[0]}_{edge[1]}", bus0=microgrid_buses.index[edge[0]], bus1=microgrid_buses.index[edge[1]])
    


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_dist import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("create_network_bis")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake) 

    n = create_network()  

    create_microgrid_network(n, snakemake.input["microgrids_buildings"])

    print(n)
    n.export_to_netcdf(snakemake.output[0])
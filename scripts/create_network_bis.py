# -*- coding: utf-8 -*-
import json

import matplotlib as plt
import pypsa
from scipy.spatial import Delaunay
from shapely.geometry import shape

# Load the GeoJSON file
with open("resources/buildings/microgrids_buildings.geojson") as f:
    data = json.load(f)

# Create a PyPSA network
network = pypsa.Network()

# Iterate over each feature in the GeoJSON file
for feature in data["features"]:
    # Get the point geometry
    point_geom = shape(feature["geometry"])

    # Create a bus at the point location with microgrid ID included in bus name
    bus_name = f"{feature['properties']['microgrid_id']}_bus_{feature['id']}"
    network.add("Bus", bus_name, x=point_geom.x, y=point_geom.y, v_nom=0.220)

# Select all buses that are part of "microgrid_1"
microgrid_1_buses = network.buses.index.str.contains("microgrid_1")
# Print the selected buses
print(microgrid_1_buses)

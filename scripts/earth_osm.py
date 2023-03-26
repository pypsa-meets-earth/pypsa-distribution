import json
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import json
import numpy as np
from scipy.spatial import Delaunay
import os

from _helpers_dist import configure_logging, sets_path_to_root


def build_plot_from_points(input_file):

    with open(input_file) as f:
        data = json.load(f)

# extract point coordinates from GeoJSON data
    points = []
    for feature in data['Data']['Node'].values():
        coords = feature['lonlat']
        points.append(coords)

# convert to numpy array
    points = np.array(points)

# create Delaunay triangulation
    tri = Delaunay(points)

    # Plot the points
    plt.plot(points[:, 0], points[:, 1], 'ko')

    # Plot the triangulation
    plt.triplot(points[:, 0], points[:, 1], tri.simplices)

    # Show the plot
    plt.show()

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_dist import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("earth_osm")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)

    configure_logging(snakemake)

    build_plot_from_points(snakemake.input["building_json"])
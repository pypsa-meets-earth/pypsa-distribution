# -*- coding: utf-8 -*-
"""
Estimates the population and the electric load of each microgrid.

Relevant Settings
-----------------
.. code:: yaml

    microgrids_list:
        microgridX: 
          lon_min:
          lon_max: 
          lat_min: 
          lat_max: 

    load:
        scaling_factor:

Inputs
------
- ``data/sample_profile.csv``: a load profile, which will be scaled through a scaling_factor to obtain the per person load

Outputs
-------
- ``resources/shapes/microgrid_shapes.geojson: a geojson file of the shape of each microgrid,
- ``resources/masked_files/masked_file_{i+1}.tif,
- ``resources/demand/microgrid_load_{i+1}.csv: the electric load of the microgid,
-------

Description
-----------

The rule :mod:`build_demand` contains functions that are used to create a shape file of the microgrid, to mask a raster with the shape file and to estimate 
the population. Then the population is multiplied for the per person load and the microgrid load is then obtained. The process applies to all the microgrids specified in config.yaml.

"""

import logging
import os

import geopandas as gpd
import pandas as pd
import rasterio
import rasterio.mask
from shapely.geometry import Polygon
from _helpers import configure_logging, get_country, sets_path_to_root
import json

_logger = logging.getLogger(__name__)
_logger.setLevel(logging.INFO)


def create_microgrid_shapes(microgrids_list, output_path): 

    """
    This function creates a rectangular shape of the microgrid and saves it as a .geojson file.
    The shape is defined by the coordinates of the angles of the rectangle.
    The resulting file is saved to the specified output_path.
    """
    
    microgrids_list=microgrids_list
    microgrids_list_df=pd.DataFrame(microgrids_list)

    microgrid_shapes = []
    microgrid_names = []

    for col in range(len(microgrids_list_df.columns)):

        values=microgrids_list_df.iloc[:, col]

        x1 = values[1] #lon_max
        y1 = values[3] #lat_max

        x2 = values[0] #lon_min
        y2 = values[3] #lat_max

        x3 = values[0] #lon_min
        y3 = values[2] #lat_min 

        x4 = values[1] #lon_max
        y4 = values[2] #lat_min


        microgrid_shape = Polygon([(x1, y1), (x2, y2), (x3, y3), (x4, y4), (x1, y1)])

        microgrid_name=f"microgrid_{col+1}"
        microgrid_shapes.append(microgrid_shape)
        microgrid_names.append(microgrid_name)

    microgrid_gdf = gpd.GeoDataFrame({'name': microgrid_names, 'geometry': microgrid_shapes})

    output_dict = json.loads(microgrid_gdf.to_json())
    output_json = json.dumps(output_dict, indent=4)

    with open(output_path, "w") as f:
        f.write(output_json)


def create_masked_file(raster_path, shapes_path, output_prefix):

    """
    Masks a raster with shapes contained in the GeoJSON file "resources/shapes/microgrid_shapes.geojson" and saves the resulting masked rasters.
    
    Parameters:
    -----------
    raster_path: str
        Path to the raster file to mask.
    shapes_path: str
        Path to the GeoJSON file containing the shapes to use for masking the raster.
    output_prefix: str
        Prefix to use for the output file names. The output files will be named 
        "{output_prefix}_{shape_index}.tif" 
    """

    # Load the GeoJSON file with the shapes to mask the raster
    shapes = gpd.read_file(shapes_path)

# Mask the raster with each shape and save each masked raster as a new file
    for i, shape in shapes.iterrows():
        with rasterio.open(raster_path) as src:
            # Mask the raster with the current shape
            masked, out_transform = rasterio.mask.mask(src, [shape.geometry], crop=True)
            out_meta = src.meta.copy()
            out_meta.update({"driver": "GTiff", "height": masked.shape[1], "width": masked.shape[2], "transform": out_transform})

        out_raster_path = f"{output_prefix}_{i+1}.tif"

        # Write the masked raster to a file
        with rasterio.open(out_raster_path, "w", **out_meta) as dest:
            dest.write(masked)   


def estimate_microgrid_population(p, sample_profile, output_files):

    """
    This function estimates the population of the microgrids based on mask files and a sample profile of electricity demand.
    """

    # Read the sample profile of electricity demand
    total_load = pd.read_csv(sample_profile)
    total_load = total_load["0"]

    # Calculate the per-person electricity demand and convert it as a pandas dataframe
    per_person_load = total_load * (1 / p)
    per_person_load = pd.DataFrame(per_person_load)

    number_microgrids = len(os.listdir("resources/masked_files"))

    for i in range(number_microgrids):
        with rasterio.open(f"resources/masked_files/country_masked_{i+1}.tif") as fp:
            data = fp.read(1)
            pop_microgrid= data[data >= 0].sum()
            
            microgrid_load  = per_person_load * pop_microgrid
            
            # Save the microgrid load to the specified output file
            output_file = os.path.join(output_files, f"microgrid_load_{i+1}.csv")
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            microgrid_load.to_csv(output_file, index=False)
   

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("build_demand")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)
    
    WorldPop=snakemake.input["WorldPop"]
    sample_profile = snakemake.input["sample_profile"]

    create_microgrid_shapes(
        snakemake.config["microgrids_list"], 
        snakemake.output["microgrid_shapes"],)

    create_masked_file(
        WorldPop, 
        snakemake.output["microgrid_shapes"],
        snakemake.output["country_masked"],
    )

    estimate_microgrid_population(
        snakemake.config["load"]["scaling_factor"],
        sample_profile,
        snakemake.output["electric_load"],
    )

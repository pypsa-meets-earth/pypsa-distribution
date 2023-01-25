"""
Estimates the population and the electric load of the microgrid.

Relevant Settings
-----------------
.. code:: yaml

    microgrids_list:
        Location:
            Centre:
                lon:
                lat:
            Sides:
                Deltalon:
                Deltalat:
        micA: 
            name:

    load:
        scaling_factor:

Inputs
------
- ``data/sample_profile.csv``: a load profile, which will be scaled through a scaling_factor to obtain the per person load

Outputs
-------
- ``data/Worldpop/population_file.tif: a tif file of the population of the selected country,
- ``resources/shapes/microgrid_shape.geojson: a geojson file of the shape of the microgrid,
- ``resources/file_dir/country_masked.tif,
- ``resources/demand/microgrid_load.csv: the electric load of the microgid,
-------

Description
-----------

The rule :mod:`build_demand` contains functions that are used to create a shape file of the microgrid, to mask a raster with the shape file and to estimate 
the population. Then the population is multiplied for the per person load and the microgrid load is then obtained.

""" 

import os
from _helpers import configure_logging, sets_path_to_root


import geopandas as gpd
import pandas as pd
import fiona
import rasterio
import rasterio.mask
import georasters as gr
import pandas as pd
import requests
import shutil
import logging
import glob


from _helpers import (
    configure_logging,
    get_country,
)

_logger = logging.getLogger(__name__)
_logger.setLevel(logging.INFO)



def create_microgrid_shape(xcenter, ycenter, DeltaX, DeltaY, name, output_path):

    """
    This function creates a rectangular shape and saves it as a .geojson file.
    The shape is defined by its center coordinates (xcenter, ycenter) and its dimensions (DeltaX, DeltaY).
    The shape is also given a name, which is saved as a property in the .geojson file.
    The resulting file is saved to the specified output_path.
    """

    # Define the coordinates of the corners of the rectangle
    x1 = xcenter - DeltaX * 0.5
    y1 = ycenter + DeltaY * 0.5

    x2 = xcenter + DeltaX * 0.5
    y2 = ycenter + DeltaY * 0.5

    x3 = xcenter - DeltaX * 0.5
    y3 = ycenter - DeltaY * 0.5

    x4 = xcenter + DeltaX * 0.5
    y4 = ycenter - DeltaY * 0.5

    microgrid_name = name
    

    my_feature = {
    "type": "FeatureCollection",
    "features": [
        {
            "type": "Feature",
            "properties": {"name": "microgrid_name"},
            "geometry": {"type": "Polygon", "coordinates": [[[x1, y1], [x2, y2], [x3, y3], [x4, y4]]]
        },
    },
    
    ],
}

    
    #Convert the feature to a GeoDataFrame
    gdf = gpd.GeoDataFrame.from_features(my_feature)

    #Save the GeoDataFrame to a .geojson file
    gdf.to_file(output_path)


def download_WorldPop_standard(
    output_path,
    country_code,
    year,
    update,
    out_logging,
    size_min,
    ):
    """
    Download tiff file for each country code using the standard method from worldpop datastore with 1kmx1km resolution.

    Parameters
    ----------
    country_code : str
        Two letter country codes of the downloaded files.
        Files downloaded from https://data.worldpop.org/ datasets WorldPop UN adjusted
    year : int
        Year of the data to download
    update : bool
        Update = true, forces re-download of files
    size_min : int
        Minimum size of each file to download
    Returns
    -------
    WorldPop_inputfile : str
        Path of the file
    WorldPop_filename : str
        Name of the file
    """

    if out_logging:
        _logger.info("Download WorldPop datasets")
    
    three_digits_code=get_country('alpha_3', alpha_2 = "SL")

    def convert_to_lowercase(string):
        return string.lower()

    three_digits_code_lower=convert_to_lowercase(three_digits_code)

    WorldPop_filename = f"{three_digits_code_lower}_ppp_{year}_constrained.tif"

    if {year} == "2019":

        WorldPop_urls = [
               f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/{year}/{three_digits_code}/{WorldPop_filename}"
        ]

    # if {year} == "2020": 

    #     WorldPop_urls = [
    #            f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/{year}/BSGM/{WorldPop_filename}"
    #     ]
    
    WorldPop_inputfile = os.path.join(
        os.getcwd(), "data", "WorldPop", f"{three_digits_code_lower}_ppp_{year}_constrained.tif" 
    ) # Input filepath tif
    

    old_file_paths = glob.glob(os.path.join("data/Worldpop", "*_ppp_*_constrained.tif"))

    for old_file_path in old_file_paths:
       os.rename(old_file_path, output_path)

    #os.remove(WorldPop_inputfile)

    if not os.path.exists(WorldPop_inputfile) or update is True:
        if out_logging:
            _logger.warning(
                f" {WorldPop_filename} does not exist, downloading to {WorldPop_inputfile}"
            )
        #  create data/osm directory
        os.makedirs(os.path.dirname(WorldPop_inputfile), exist_ok=True)

        loaded = False

        WorldPop_urls = [
        f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/{year}/{three_digits_code}/{WorldPop_filename}",
        ]
        #f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/{year}/{three_digits_code}/{WorldPop_filename}"

    
        
        for WorldPop_url in WorldPop_urls:
            with requests.get(WorldPop_url, stream=True) as r:
                with open(WorldPop_inputfile, "wb") as f:
                    if float(r.headers["Content-length"]) > size_min:
                        shutil.copyfileobj(r.raw, f)
                        loaded = True
                        break
        if not loaded:
            _logger.error(f" Impossible to download {WorldPop_filename}")



def create_masked_file(raster_path, geojson_path, output_path):

    """
    This function masks a raster with a shape defined in a geojson file.
    The raster file is specified with the raster_path, the shape is specified with the geojson_path,
    and the resulting masked raster is saved to the specified output_path.
    """
    
    # Read the geojson file and convert it to a shapefile
    gdf = gpd.read_file(geojson_path)
    gdf.to_file('microgrid_shape.shp')

    # Open the shapefile and extract the shape geometry
    with fiona.open("microgrid_shape.shp", "r") as shapefile:
      shapes = [feature["geometry"] for feature in shapefile]


    # Open the raster and mask it using the shapes
    with rasterio.open(raster_path) as src:
        out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
        out_meta = src.meta

    # update the metadata for the output raster
    out_meta.update({"driver": "GTiff",
                     "height": out_image.shape[1],
                     "width": out_image.shape[2],
                     "transform": out_transform})

    # Save the masked raster to the specified output path
    with rasterio.open(output_path,"w", **out_meta) as dest:
        dest.write(out_image)


 
def estimate_microgrid_population(masked_file, p, sample_profile, output_file):

    """
    This function estimates the population of a microgrid based on a mask file and a sample profile of electricity demand.
    The mask file is specified with the masked_file, the sample profile is specified with the sample_profile,
    and the output file for the estimated microgrid population is specified with the output_file.
    """
   
    pop_microgrid = gr.from_file(masked_file)
    
    #Read the pop_mircrogid file and convert it to a geopandas dataframe
    pop_microgrid=pop_microgrid.to_geopandas() 

    # Sum the values of the mask file to get the total population of the microgrid
    pop_microgrid=(pop_microgrid['value'].sum())  

    # Read the sample profile of electricity demand
    total_load=pd.read_csv(sample_profile)
    total_load = total_load["0"]

    # Calculate the per-person electricity demand and convert it as a pandas dataframe
    per_person_load=total_load*(1/p) 
    per_person_load=pd.DataFrame(per_person_load)

    #Calculate the microgrid electric load
    microgrid_load=per_person_load*pop_microgrid 
    
    # Save the microgrid load to the specified output file
    microgrid_load.to_csv(output_file, index=False)
    
    return microgrid_load
    


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("build_demand")
        configure_logging(snakemake)

    sample_profile=snakemake.input["sample_profile"]

    create_microgrid_shape(
        snakemake.config["microgrids_list"]["Location"]["Centre"]["lon"],
        snakemake.config["microgrids_list"]["Location"]["Centre"]["lat"],
        snakemake.config["microgrids_list"]["Location"]["Sides"]["Deltalon"],
        snakemake.config["microgrids_list"]["Location"]["Sides"]["Deltalat"],
        snakemake.config["microgrids_list"]["micA"]["name"],
        snakemake.output["microgrid_shape"],
    )

    download_WorldPop_standard(
    snakemake.output["Worldpop_data"],
    snakemake.config["countries"],
    snakemake.config["year"],
    False,
    False,
    300,
)

    create_masked_file(f"data/Worldpop/population_file.tif",
                    f"resources/shapes/microgrid_shape.geojson",
                    snakemake.output["country_masked"])

    estimate_microgrid_population(f"resources/file_dir/country_masked.tif", 
                                snakemake.config["load"]["scaling_factor"],       
                                sample_profile, 
                                snakemake.output["electric_load"])

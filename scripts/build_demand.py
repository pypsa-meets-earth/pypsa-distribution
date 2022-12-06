import os
from _helpers import configure_logging, sets_path_to_root


import geopandas as gpd
import pandas as pd
import fiona
import rasterio
import rasterio.mask
import georasters as gr
import pandas as pd
import geojson
import shutil
#%%
# #This script downloads WorldPop data for SL country. 2019 data are taken, since 2020 data are not available for Sierra Leone

# import requests 
# import shutil
# import os

# def download_WorldPop_standard(
#     country_code,
#     year=2019,
#     update=False,
#     size_min=300
# ):
#     """
#     Download tiff file for each country code using the standard method from worldpop datastore with 1kmx1km resolution.

#     Parameters
#     ----------
#     country_code : str
#         Two letter country codes of the downloaded files.
        # Files downloaded from https://data.worldpop.org/ datasets WorldPop UN adjusted
#     year : int
#         Year of the data to download
#     update : bool
#         Update = true, forces re-download of files
#     size_min : int
#         Minimum size of each file to download
#     Returns
#     -------
#     WorldPop_inputfile : str
#         Path of the file
#     WorldPop_filename : str
#         Name of the file
#     """
#     WorldPop_filename = f"sle_ppp_{year}_constrained.tif"
#     WorldPop_urls = [
#             f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2019/BSGM/SLE/{WorldPop_filename}"]
#     WorldPop_inputfile = os.path.join(
#         os.getcwd(), "data", "WorldPop", WorldPop_filename
#     )  # Input filepath tif

#     if not os.path.exists(WorldPop_inputfile) or update is True:
    
#         #  create data/osm directory
#         os.makedirs(os.path.dirname(WorldPop_inputfile), exist_ok=True)

#         loaded = False
#         for WorldPop_url in WorldPop_urls:
#             with requests.get(WorldPop_url, stream=True) as r:
#                 with open(WorldPop_inputfile, "wb") as f:
#                     if float(r.headers["Content-length"]) > size_min:
#                         shutil.copyfileobj(r.raw, f)
#                         loaded = True
#                         break

#     return WorldPop_inputfile, WorldPop_filename


# download_WorldPop_standard(
#                     config["countries"], config["year"], False, 300)


#****

#The vertexes of the rectangular are defined: they are the point P1=(x1,y1), P2=(x2,y2), P3=(x3,y3), P4=(x4,y4)
#x stands for longitude, y stands for latitude

def create_microgrid_shape(xcenter, ycenter, DeltaX, DeltaY, name):

    x1 = xcenter - DeltaX*0.5
    y1 = ycenter + DeltaY*0.5

    x2 = xcenter + DeltaX*0.5
    y2 = ycenter + DeltaY*0.5

    x3 = xcenter - DeltaX*0.5 
    y3 = ycenter - DeltaY*0.5

    x4 = xcenter + DeltaX*0.5
    y4 = ycenter - DeltaY*0.5

    microgrid_name=name
    my_feature=[]

    my_feature={
      "type": "Feature",
      "geometry": {
        "type": "Polygon",
        "coordinates":  [
            [
              [x1,y1],
              [x2,y2],
              [x3,y3],         
              [x4,y4],
            ]
        ]
      },

      "properties": {
        "Microgrid": [microgrid_name]
      }
    }

    return(my_feature)

# my_feature is converted into a .geojson file

def writeToGeojsonFile(path, fileName, data):
    filePathNameWExt = './' + path + '/' + fileName + '.geojson'
    with open(filePathNameWExt, 'w') as fp:
        geojson.dump(data, fp)
  
    geojson_filename=f"microgrid_shape.geojson"

    source = os.path.join(
        os.getcwd(), geojson_filename)
    destination = os.path.join(
        os.path.abspath(os.curdir), "resources", "shapes"
    )  

    shutil.copy(source, destination)



def from_geojson_to_tif():

    gdf = gpd.read_file('microgrid_shape.geojson')
    gdf.to_file('microgrid_shape.shp')

    with fiona.open("microgrid_shape.shp", "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]
    
    with rasterio.open(f"data/Worldpop/sle_ppp_2019_constrained.tif") as src:
        out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
        out_meta = src.meta

    out_meta.update({"driver": "GTiff",
                     "height": out_image.shape[1],
                     "width": out_image.shape[2],
                    "transform": out_transform})

    with rasterio.open("SL.masked.tif", "w", **out_meta) as dest:
        dest.write(out_image)

def estimate_microgrid_population():
    myRaster = f"data/Worldpop/sle_ppp_2019_constrained.tif"
    total_pop= gr.from_file(myRaster)
    
    total_pop=total_pop.to_geopandas() 

    total_pop=(total_pop['value'].sum()) #Total SL population

    myRaster = 'SL.masked.tif'
    pop_microgrid = gr.from_file(myRaster)
    
    pop_microgrid=pop_microgrid.to_geopandas() 

    pop_microgrid=(pop_microgrid['value'].sum()) #Microgrid population

#I import the dataframe of electricity demand for Africa
    df_demand=pd.read_excel(f"data/Africa.xlsx", index_col = None)

#I select the rows related to Benin (since there are no data for SL) 
    df_demand_SL=df_demand.loc[26280:35039, :]

# I select the column "electricity demand"
    df_demand_SL=df_demand_SL["Electricity demand"]

    df_demand_SL=pd.DataFrame(df_demand_SL)

    p=(pop_microgrid/total_pop)*100 #Coefficient

    electric_load=df_demand_SL/p #Electric load of the minigrid

    electric_load=electric_load.to_excel('electric_load.xlsx', index=False) 
    
    xlsx_filename=f"electric_load.xlsx"

    source = os.path.join(
        os.getcwd(), xlsx_filename)
    destination = os.path.join(
        os.path.abspath(os.curdir), "resources", "demand"
    )  #TODO: demand folder that creates automatically

    shutil.copy(source, destination)

    return electric_load


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("build_demand")

    configure_logging(snakemake)

    out = snakemake.output

    my_feature=create_microgrid_shape(
        snakemake.config["microgrids_list"]["Location"]["Centre"]["lon"],
        snakemake.config["microgrids_list"]["Location"]["Centre"]["lat"],
        snakemake.config["microgrids_list"]["Location"]["Sides"]["Deltalon"],
        snakemake.config["microgrids_list"]["Location"]["Sides"]["Deltalat"],
        snakemake.config["microgrids_list"]["micA"]["name"],
    )

    writeToGeojsonFile('./','microgrid_shape', my_feature)

    from_geojson_to_tif()
    
    electric_load=estimate_microgrid_population()
    #%%
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

    if not os.path.exists(os.path.join(os.getcwd(), "resources", "shapes")):
        os.makedirs("resources/shapes")
   
    filePathNameWExt = './' + path + '/resources/shapes/' + fileName + '.geojson'
    with open(filePathNameWExt, 'w') as fp:
        geojson.dump(data, fp)


def from_geojson_to_tif():

    gdf = gpd.read_file(f"resources/shapes/microgrid_shape.geojson")
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

    demand_microgrid=df_demand_SL/p #Electric load of the minigrid

    return demand_microgrid
    


def create_load_file():
    
    if not os.path.exists(os.path.join(os.getcwd(), "resources", "demand")):
        os.makedirs("resources/demand")

    electric_load=demand_microgrid.to_excel('electric_load.xlsx', index=False) 

    xlsx_filename=f"electric_load.xlsx"

      
    shutil.copy(os.path.join(
        os.getcwd(), xlsx_filename), os.path.join(
        os.path.abspath(os.curdir), "resources", "demand") )
        
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
    
    demand_microgrid=estimate_microgrid_population()

    electric_load_file=create_load_file()
    #%%
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
import pypsa


def create_microgrid_shape(xcenter, ycenter, DeltaX, DeltaY, name, output_path):
    # define the coordinates of the rectangle
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
            "properties": {"col1": "name1"},
            "geometry": {"type": "Polygon", "coordinates": [[[x1, y1], [x2, y2], [x3, y3], [x4, y4]]]
        },
    },
    
    ],
}
    
    #my_feature is converted into a .geojson file 
    gdf = gpd.GeoDataFrame.from_features(my_feature)
    gdf.to_file(output_path)


def create_masked_file(WorldPop_data):

    if not os.path.exists(os.path.join(os.getcwd(), "resources", "file_dir")):
        os.makedirs("resources/file_dir") #I created a directory for the following files to be stored

    gdf = gpd.read_file(f"resources/shapes/microgrid_shape.geojson")
    gdf.to_file('./resources/file_dir/microgrid_shape.shp')

    with fiona.open(f"resources/file_dir/microgrid_shape.shp", "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]
    
    with rasterio.open(WorldPop_data) as src:
        out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
        out_meta = src.meta

    out_meta.update({"driver": "GTiff",
                     "height": out_image.shape[1],
                     "width": out_image.shape[2],
                    "transform": out_transform})
   
    with rasterio.open("./resources/file_dir/SL.masked.tif", "w", **out_meta) as dest:
        dest.write(out_image)

    

def estimate_microgrid_population(sample_profile):

    pop_microgrid = gr.from_file(f"resources/file_dir/SL.masked.tif")
    
    pop_microgrid=pop_microgrid.to_geopandas() 

    pop_microgrid=(pop_microgrid['value'].sum())  #Microgrid population

    #I import the sample_profile file
    total_load=pd.read_csv(sample_profile)
    total_load=total_load["0"]
    per_person_load=total_load*(1/1500) #file of electricity demand per-person
    
    per_person_load=pd.DataFrame(per_person_load)
    microgrid_load=per_person_load*pop_microgrid #Electric load of the microgrid
        
    if not os.path.exists(os.path.join(os.getcwd(), "resources", "demand")):
        os.makedirs("resources/demand")

    microgrid_load=microgrid_load.to_csv("./resources/demand/microgrid_load.csv", index=False) 
    
    return microgrid_load


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("build_demand")
        configure_logging(snakemake)


    # WorldPop_data=snakemake.input["WorldPop"]
    # sample_profile=snakemake.input["sample_profile"]

    create_microgrid_shape(
        snakemake.config["microgrids_list"]["Location"]["Centre"]["lon"],
        snakemake.config["microgrids_list"]["Location"]["Centre"]["lat"],
        snakemake.config["microgrids_list"]["Location"]["Sides"]["Deltalon"],
        snakemake.config["microgrids_list"]["Location"]["Sides"]["Deltalat"],
        snakemake.config["microgrids_list"]["micA"]["name"],
        snakemake.output["microgrid_shape"],
    )

    

    # create_masked_file(WorldPop_data)

    # microgrid_load=estimate_microgrid_population(sample_profile)
 

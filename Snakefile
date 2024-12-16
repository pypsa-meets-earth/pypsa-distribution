import sys

sys.path.append("./scripts")
sys.path.append("./pypsa-earth/scripts")

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.io import expand
from build_test_configs import create_test_config
from _helpers import create_country_list
from os.path import normpath, exists, isdir
from shutil import copyfile
from pathlib import Path
from _helpers import create_country_list

import sys

sys.path.append("./scripts")

from os.path import normpath, exists, isdir

HTTP = HTTPRemoteProvider()


COSTS = "data/costs.csv"
PROFILE = "data/sample_profile.csv"

PYPSAEARTH_FOLDER = "pypsa-earth"

if "config" not in globals() or not config:  # skip when used as sub-workflow
    if not exists("config.yaml"):
        # # prepare pypsa-earth config
        # create_test_config(
        #     "./config.pypsa-earth.yaml", "./config.distribution.yaml", "./config.yaml"
        # )
        copyfile("config.distribution.yaml", "config.yaml")

    configfile: "config.pypsa-earth.yaml"
    configfile: "config.yaml"


config["countries"] = create_country_list(config["countries"])


config["countries"] = create_country_list(config["countries"])

run = config.get("run", {})
RDIR = run["name"] + "/" if run.get("name") else ""
countries = config["countries"]

ATLITE_NPROCESSES = config["atlite"].get("nprocesses", 5)


wildcard_constraints:
    ll="[a-z0-9\.]+",
    simpl="[a-zA-Z0-9]*|all",
    clusters="[0-9]+m?|all",
    opts="[-+a-zA-Z0-9]*",
    sopts="[-+a-zA-Z0-9\.\s]*",
    user_type="[a-zA-Z0-9]*",


if not config.get("disable_subworkflow", False):

    subworkflow pypsaearth:
        workdir:
            PYPSAEARTH_FOLDER
        snakefile:
            PYPSAEARTH_FOLDER + "/Snakefile"
        configfile:
            "./config.pypsa-earth.yaml"


if config.get("disable_subworkflow", False):

    def pypsaearth(path):
        return PYPSAEARTH_FOLDER + "/" + path


# rule clean:
#     run:
#         shell("snakemake -j 1 solve_network --delete-all-output")


rule ramp_build_demand_profile:
    params:
        ramp=config["ramp"],
        snapshoots=config["snapshots"],
    input:
        user_description="data/ramp/{user_type}.xlsx",
    output:
        daily_demand_profiles="resources/ramp/daily_demand_{user_type}.xlsx",
        daily_type_demand_profile="resources/ramp/daily_type_demand_{user_type}.xlsx",
    log:
        "logs/ramp_build_demand_profile_{user_type}.log",
    benchmark:
        "benchmarks/ramp_build_demand_profile_{user_type}"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "scripts/ramp_build_demand_profile.py"


rule build_demand:
    params:
        tier=config["tier"],
        snapshots=config["snapshots"],
        build_demand_model=config["build_demand_type"],
    input:
        **{
            f"profile_{user_file.stem}": f"resources/ramp/daily_type_demand_{user_file.stem}.xlsx"
            for user_file in Path("data/ramp/").glob("[a-zA-Z0-9]*.xlsx")
        },
        sample_profile=PROFILE,
        building_csv="resources/buildings/buildings_type.csv",
        create_network="networks/base.nc",
        microgrid_shapes="resources/shapes/microgrid_shapes.geojson",
        clusters_with_buildings="resources/buildings/cluster_with_buildings.geojson",
    output:
        electric_load="resources/demand/microgrid_load.csv",
    log:
        "logs/build_demand.log",
    benchmark:
        "benchmarks/build_demand"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "scripts/build_demand.py"


rule build_shapes:
    output:
        microgrid_shapes="resources/shapes/microgrid_shapes.geojson",
        microgrid_bus_shapes="resources/shapes/microgrid_bus_shapes.geojson",
    log:
        "logs/build_shapes.log",
    benchmark:
        "benchmarks/build_shapes"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "scripts/build_shapes.py"


rule create_network:
    input:
        clusters="resources/buildings/clustered_buildings.geojson",
    output:
        "networks/base.nc",
    log:
        "logs/create_network.log",
    benchmark:
        "benchmarks/create_network"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "scripts/create_network.py"


if config["enable"].get("download_osm_buildings", True):

    rule download_osm_data:
        output:
            building_resources="resources/" + RDIR + "osm/raw/all_raw_buildings.geojson",
        log:
            "logs/" + RDIR + "download_osm_data.log",
        benchmark:
            "benchmarks/" + RDIR + "download_osm_data"
        threads: 1
        resources:
            mem_mb=3000,
        script:
            "scripts/download_osm_data.py"


rule clean_earth_osm_data:
    input:
        all_buildings="resources/" + RDIR + "osm/raw/all_raw_buildings.geojson",
        microgrid_shapes="resources/shapes/microgrid_shapes.geojson",
    output:
        microgrid_building="resources/buildings/microgrid_building.geojson",
    log:
        "logs/clean_earth_osm_data.log",
    benchmark:
        "benchmarks/clean_earth_osm_data"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "scripts/clean_earth_osm_data.py"


rule cluster_buildings:
    params:
        crs=config["crs"],
        house_area_limit=config["house_area_limit"],
    input:
        buildings_geojson="resources/buildings/microgrid_building.geojson",
    output:
        clusters="resources/buildings/clustered_buildings.geojson",
        clusters_with_buildings="resources/buildings/cluster_with_buildings.geojson",
        buildings_type="resources/buildings/buildings_type.csv",
    log:
        "logs/cluster_buildings.log",
    benchmark:
        "benchmarks/cluster_buildings"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "scripts/cluster_buildings.py"


rule build_renewable_profiles:
    params:
        crs=config["crs"],
        renewable=config["renewable"],
        countries=config["countries"],
        alternative_clustering=config["cluster_options"]["alternative_clustering"],
    input:
        natura=pypsaearth("resources/natura.tiff"),
        copernicus=pypsaearth(
            "data/copernicus/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif"
        ),
        gebco=pypsaearth("data/gebco/GEBCO_2021_TID.nc"),
        country_shapes="resources/shapes/microgrid_shapes.geojson",
        offshore_shapes=pypsaearth("resources/shapes/offshore_shapes.geojson"),
        hydro_capacities="pypsa-earth/data/hydro_capacities.csv",
        eia_hydro_generation="pypsa-earth/data/eia_hydro_annual_generation.csv",
        powerplants="resources/powerplants.csv",
        regions="resources/shapes/microgrid_bus_shapes.geojson",
        cutout=lambda w: pypsaearth(
            "cutouts/" + config["renewable"][w.technology]["cutout"] + ".nc"
        ),
    output:
        profile="resources/renewable_profiles/profile_{technology}.nc",
    log:
        "logs/build_renewable_profile_{technology}.log",
    benchmark:
        "benchmarks/build_renewable_profiles_{technology}"
    threads: ATLITE_NPROCESSES
    resources:
        mem_mb=ATLITE_NPROCESSES * 5000,
    script:
        pypsaearth("scripts/build_renewable_profiles.py")


rule add_electricity:
    input:
        **{
            f"profile_{tech}": f"resources/renewable_profiles/profile_{tech}.nc"
            for tech in config["tech_modelling"]["general_vre"]
        },
        create_network="networks/base.nc",
        tech_costs=COSTS,
        load_file="resources/demand/microgrid_load.csv",
        powerplants="resources/powerplants.csv",
    output:
        "networks/elec.nc",
    log:
        "logs/add_electricity.log",
    benchmark:
        "benchmarks/add_electricity"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "scripts/add_electricity.py"


# if config["monte_carlo"]["options"].get("add_to_snakefile", False) == False:

# rule solve_network:
#     input:
#         "networks/elec.nc",
#     output:
#         "networks/results/elec.nc",
#     log:
#         "logs/solve_network.log",
#     benchmark:
#         "benchmarks/solve_network"
#     threads: 1
#     resources:
#         mem_mb=3000,
#     script:
#         "scripts/solve_network.py"


rule solve_network:
    input:
        "networks/elec.nc",
    output:
        "networks/results/elec.nc",
    log:
        "logs/solve_network.log",
    benchmark:
        "benchmarks/solve_network"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "scripts/solve_network.py"

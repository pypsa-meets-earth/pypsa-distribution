import sys

sys.path.append("./scripts")
sys.path.append("./pypsa-earth/scripts")

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.io import expand
from build_test_configs import create_test_config
from _helpers import create_country_list
from os.path import normpath, exists, isdir
from shutil import copyfile

import sys

sys.path.append("./scripts")

from os.path import normpath, exists, isdir

HTTP = HTTPRemoteProvider()


COSTS = "data/costs.csv"
PROFILE = "data/sample_profile.csv"

if "config" not in globals() or not config:  # skip when used as sub-workflow
    if not exists("config.yaml"):
        # prepare pypsa-earth config
        create_test_config(
            "./config.pypsa-earth.yaml", "./config.distribution.yaml", "./config.yaml"
        )
        # copyfile("config.distribution.yaml", "config.yaml")

    configfile: "config.yaml"

# convert country list according to the desired region
config["countries"] = create_country_list(config["countries"])

run = config.get("run", {})
RDIR = run["name"] + "/" if run.get("name") else ""

ATLITE_NPROCESSES = config["atlite"].get("nprocesses", 5)


wildcard_constraints:
    ll="[a-z0-9\.]+",
    simpl="[a-zA-Z0-9]*|all",
    clusters="[0-9]+m?|all",
    opts="[-+a-zA-Z0-9]*",
    sopts="[-+a-zA-Z0-9\.\s]*",
    discountrate="[-+a-zA-Z0-9\.\s]*",


subworkflow pypsaearth:
    workdir:
        "./pypsa-earth"
    snakefile:
        "./pypsa-earth/Snakefile"
    configfile:
        "./config.yaml"


rule clean:
    run:
        shell("snakemake -j 1 solve_network --delete-all-output")


rule build_demand:
    input:
        sample_profile=PROFILE,
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


rule clean_earth_osm_data:
    input:
        #buildings_json="resources/buildings/buildings.json",
        microgrid_shapes="resources/shapes/microgrid_shapes.geojson",
    output:
        buildings_geojson="resources/buildings/buildings.geojson",
        microgrids_buildings="resources/buildings/microgrids_buildings.geojson",
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
    input:
        buildings_geojson="resources/buildings/buildings.geojson",
    output:
        cleaned_buildings_geojson="resources/buildings/cleaned_buildings.geojson",
        clusters="resources/buildings/clustered_buildings.geojson",
        clusters_with_buildings="resources/buildings/cluster_with_buildings.geojson",
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
        country_shapes=pypsaearth("resources/shapes/country_shapes.geojson"),
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


if config["enable"].get("download_osm_data", True):

    rule download_osm_data:
        params:
            countries=config["countries"],
        output:
            "resources", RDIR, "osm", "raw"
        log:
            "logs/" + RDIR + "download_osm_data.log",
        benchmark:
            "benchmarks/" + RDIR + "download_osm_data"
        script:
            "scripts/download_osm_data.py"

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
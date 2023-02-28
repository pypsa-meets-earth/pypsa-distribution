import sys

sys.path.append("./scripts")

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.io import expand
from _helpers import merge_yamls
from os.path import normpath, exists, isdir

import sys

sys.path.append("./scripts")

from os.path import normpath, exists, isdir

HTTP = HTTPRemoteProvider()


configfile: "config.yaml"


COSTS = "data/costs.csv"
PROFILE = "data/sample_profile.csv"

# prepare pypsa-earth config
merge_yamls(
    "./pypsa-earth/config.default.yaml", "./config.yaml", "./config.pypsa-earth.yaml"
)

# ATLITE_NPROCESSES = config["atlite"].get("nprocesses", 20)


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


rule build_demand:
    input:
        WorldPop="data/WorldPop/population_file.tif",
        sample_profile=PROFILE,
    output:
        microgrid_shapes="resources/shapes/microgrid_shapes.geojson",
        country_masked="resources/masked_files/country_masked",
        electric_load="resources/demand",
    log:
        "logs/build_demand.log",
    benchmark:
        "benchmarks/build_demand"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "scripts/build_demand.py"

rule create_network:
    output:
        # [f"networks/base_{i+1}.nc" for i in range(2)]
        [f"networks/base_{i+1}.nc" for i in range(len(config['microgrids_list']))]
    log:
        "logs/create_network.log",
    benchmark:
        "benchmarks/create_network"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "scripts/create_network.py"

# rule create_network:
#     output:
#         #"networks/base.nc",
#         #[f"networks/base_{i}.nc" for i in range(3)],
#         [f"networks/base_{i+1}.nc" for i in range(0,2)],

#     log:
#         "logs/create_network.log",
#     benchmark:
#         "benchmarks/create_network"
#     threads: 1
#     resources:
#         mem_mb=3000,
#     script:
#         "scripts/create_network.py"


# rule build_renewable_profiles:
#     input:
#         base_network="networks/base.nc",
#         natura=pypsaearth("resources/natura.tiff"),
#         copernicus=pypsaearth(
#             "data/copernicus/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif"
#         ),
#         gebco=pypsaearth("data/gebco/GEBCO_2021_TID.nc"),
#         country_shapes="resources/shapes/country_shapes.geojson",
#         offshore_shapes="resources/shapes/offshore_shapes.geojson",
#         hydro_capacities=pypsaearth("data/hydro_capacities.csv"),
#         eia_hydro_generation=pypsaearth("data/eia_hydro_annual_generation.csv"),
#         powerplants="resources/powerplants.csv",
#         regions=lambda w: (
#             "resources/bus_regions/regions_onshore.geojson"
#             if w.technology in ("onwind", "solar", "hydro")
#             else "resources/bus_regions/regions_offshore.geojson"
#         ),
#         cutout=lambda w: "cutouts/"
#         + config["renewable"][w.technology]["cutout"]
#         + ".nc",
#     output:
#         profile="resources/renewable_profiles/profile_{technology}.nc",
#     log:
#         "logs/build_renewable_profile_{technology}.log",
#     benchmark:
#         "benchmarks/build_renewable_profiles_{technology}"
#     threads: ATLITE_NPROCESSES
#     resources:
#         mem_mb=ATLITE_NPROCESSES * 5000,
#     script:
#         pypsaearth("scripts/build_renewable_profiles.py")


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

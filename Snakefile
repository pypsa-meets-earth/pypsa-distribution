import sys

sys.path.append("./scripts")

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from _helpers import merge_yamls

HTTP = HTTPRemoteProvider()


configfile: "config.yaml"


COSTS = "data/costs.csv"


# prepare pypsa-earth config
merge_yamls(
    "./pypsa-earth/config.default.yaml", "./config.yaml", "./config.pypsa-earth.yaml"
)

ATLITE_NPROCESSES = config["atlite"].get("nprocesses", 20)


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


rule build_shapes:
    output:
        "resources/shapes/shapes.geojson",
    log:
        "logs/build_shapes.log",
    benchmark:
        "benchmarks/build_shapes"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "scripts/build_shapes.py"


rule build_demand:
    input:
        shapes="resources/shapes/shapes.geojson",
    output:
        "resources/demand/electric_load.xlsx",
    log:
        "logs/build_demand.log",
    benchmark:
        "benchmarks/build_demand"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "scripts/build_demand.py"


rule base_network:
    input:
        shapes="resources/shapes/shapes.geojson",
    output:
        "networks/base.nc",
    log:
        "logs/base_network.log",
    benchmark:
        "benchmarks/base_network"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "scripts/base_network.py"


rule build_renewable_profiles:
    input:
        base_network="networks/base.nc",
        natura=pypsaearth("resources/natura.tiff"),
        copernicus=pypsaearth(
            "data/copernicus/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif"
        ),
        gebco=pypsaearth("data/gebco/GEBCO_2021_TID.nc"),
        country_shapes="resources/shapes/country_shapes.geojson",
        offshore_shapes="resources/shapes/offshore_shapes.geojson",
        hydro_capacities=pypsaearth("data/hydro_capacities.csv"),
        eia_hydro_generation=pypsaearth("data/eia_hydro_annual_generation.csv"),
        powerplants="resources/powerplants.csv",
        regions=lambda w: (
            "resources/bus_regions/regions_onshore.geojson"
            if w.technology in ("onwind", "solar", "hydro")
            else "resources/bus_regions/regions_offshore.geojson"
        ),
        cutout=lambda w: "cutouts/"
        + config["renewable"][w.technology]["cutout"]
        + ".nc",
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
        base_network="networks/base.nc"
        ** {
            f"profile_{tech}": f"resources/renewable_profiles/profile_{tech}.nc"
            for tech in config["renewable"]
            if tech in config["electricity"]["renewable_carriers"]
        },
        demand="resources/demand/electric_load.xlsx",  # path for the nodal demand
        tech_costs=COSTS,
        regions="resources/bus_regions/regions_onshore.geojson",
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

from retrieve_databundle_light import (
    datafiles_retrivedatabundle,
    get_best_bundles_in_snakemake,
)

PYPSAEARTH_FOLDER = "../pypsa-distribution/pypsa-earth"

if not config.get("disable_subworkflow", False):

    subworkflow pypsaearth:
        workdir:
            PYPSAEARTH_FOLDER
        snakefile:
            PYPSAEARTH_FOLDER + "/Snakefile"
        configfile:
            "config.pypsa-earth.yaml"


if config.get("disable_subworkflow", False):

    def pypsaearth(path):
        return PYPSAEARTH_FOLDER + "/" + path


COSTS = "data/costs.csv"
PROFILE = "data/sample_profile.csv"


configfile: "config.pypsa-earth.yaml"


configfile: "pypsa-earth/configs/bundle_config.yaml"


rule dist_ramp_build_demand_profile:
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
        "../scripts/dist_ramp_build_demand_profile.py"


rule dist_build_demand:
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
        microgrid_shapes="resources/shapes/microgrid_shapes.geojson",
        clusters_with_buildings="resources/buildings/cluster_with_buildings.geojson",
    output:
        electric_load="resources/demand/microgrid_load.csv",
    log:
        "logs/dist_build_demand.log",
    benchmark:
        "benchmarks/dist_build_demand"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "../scripts/dist_build_demand.py"


rule dist_build_shapes:
    params:
        countries=config["countries"],
    output:
        microgrid_shapes="resources/shapes/microgrid_shapes.geojson",
        microgrid_bus_shapes="resources/shapes/microgrid_bus_shapes.geojson",
    log:
        "logs/dist_build_shapes.log",
    benchmark:
        "benchmarks/dist_build_shapes"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "../scripts/dist_build_shapes.py"


if config.get("mode") != "brown_field":

    rule dist_cluster_buildings:
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
            "logs/dist_cluster_buildings.log",
        benchmark:
            "benchmarks/dist_cluster_buildings"
        threads: 1
        resources:
            mem_mb=3000,
        script:
            "../scripts/dist_cluster_buildings.py"

    rule dist_create_network:
        input:
            clusters="resources/buildings/clustered_buildings.geojson",
            load="resources/demand/microgrid_load.csv",
        output:
            "networks/" + RDIR + "base.nc",
        log:
            "logs/dist_create_network.log",
        benchmark:
            "benchmarks/dist_create_network"
        threads: 1
        resources:
            mem_mb=3000,
        script:
            "../scripts/dist_create_network.py"


if config["enable"].get("download_osm_buildings", True):

    rule dist_download_osm_data:
        output:
            buildings_resources="resources/"
            + RDIR
            + "osm/raw/all_raw_buildings.geojson",
            lines_resources="resources/" + RDIR + "osm/raw/all_raw_lines.geojson",
            cables_resources="resources/" + RDIR + "osm/raw/all_raw_cables.geojson",
            generators_resources="resources/"
            + RDIR
            + "osm/raw/all_raw_generators.geojson",
            substations_resources="resources/"
            + RDIR
            + "osm/raw/all_raw_substations.geojson",
            poles_resources="resources/" + RDIR + "osm/raw/all_raw_poles.geojson",
        log:
            "logs/" + RDIR + "dist_download_osm_data.log",
        benchmark:
            "benchmarks/" + RDIR + "dist_download_osm_data"
        threads: 1
        resources:
            mem_mb=3000,
        script:
            "../scripts/dist_download_osm_data.py"


rule dist_clean_earth_osm_data:
    input:
        all_buildings="resources/" + RDIR + "osm/raw/all_raw_buildings.geojson",
        microgrid_shapes="resources/shapes/microgrid_shapes.geojson",
    output:
        microgrid_building="resources/buildings/microgrid_building.geojson",
    log:
        "logs/dist_clean_earth_osm_data.log",
    benchmark:
        "benchmarks/dist_clean_earth_osm_data"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "../scripts/dist_clean_earth_osm_data.py"


if config.get("mode") == "brown_field":

    rule dist_clean_osm_data:
        params:
            crs=config["crs"],
            clean_osm_data_options=config["clean_osm_data_options"],
        input:
            cables="resources/" + RDIR + "osm/raw/all_raw_cables.geojson",
            generators="resources/" + RDIR + "osm/raw/all_raw_generators.geojson",
            lines="resources/" + RDIR + "osm/raw/all_raw_lines.geojson",
            substations="resources/" + RDIR + "osm/raw/all_raw_substations.geojson",
            country_shapes="resources/shapes/microgrid_shapes.geojson",
            offshore_shapes=pypsaearth("resources/shapes/offshore_shapes.geojson"),
            africa_shape=pypsaearth("resources/shapes/africa_shape.geojson"),
        output:
            generators="resources/" + RDIR + "osm/clean/all_clean_generators.geojson",
            generators_csv="resources/" + RDIR + "osm/clean/all_clean_generators.csv",
            lines="resources/" + RDIR + "osm/clean/all_clean_lines.geojson",
            substations="resources/" + RDIR + "osm/clean/all_clean_substations.geojson",
        log:
            "logs/" + RDIR + "clean_osm_data.log",
        benchmark:
            "benchmarks/" + RDIR + "clean_osm_data"
        script:
            pypsaearth("scripts/clean_osm_data.py")

    rule dist_build_osm_network:
        params:
            build_osm_network=config.get("build_osm_network", {}),
            countries=config["countries"],
            crs=config["crs"],
        input:
            generators="resources/" + RDIR + "osm/clean/all_clean_generators.geojson",
            lines="resources/" + RDIR + "osm/clean/all_clean_lines.geojson",
            substations="resources/" + RDIR + "osm/clean/all_clean_substations.geojson",
            country_shapes="resources/" + RDIR + "shapes/microgrid_shapes.geojson",
        output:
            lines="resources/" + RDIR + "base_network/all_lines_build_network.csv",
            converters="resources/"
            + RDIR
            + "base_network/all_converters_build_network.csv",
            transformers="resources/"
            + RDIR
            + "base_network/all_transformers_build_network.csv",
            substations="resources/" + RDIR + "base_network/all_buses_build_network.csv",
        log:
            "logs/" + RDIR + "dist_build_osm_network.log",
        benchmark:
            "benchmarks/" + RDIR + "dist_build_osm_network"
        script:
            "../scripts/dist_build_osm_network.py"

    rule dist_cluster_buildings:
        params:
            crs=config["crs"],
            house_area_limit=config["house_area_limit"],
            voltage_node_cluster=config["electricity"]["voltage_node_cluster"],
        input:
            buildings_geojson="resources/buildings/microgrid_building.geojson",
            all_nodes_brown_field="resources/"
            + RDIR
            + "base_network/all_buses_build_network.csv",
        output:
            clusters="resources/buildings/clustered_buildings.geojson",
            clusters_with_buildings="resources/buildings/cluster_with_buildings.geojson",
            buildings_type="resources/buildings/buildings_type.csv",
        log:
            "logs/dist_cluster_buildings.log",
        benchmark:
            "benchmarks/dist_cluster_buildings"
        threads: 1
        resources:
            mem_mb=3000,
        script:
            "../scripts/dist_cluster_buildings.py"

    rule dist_base_network:
        params:
            voltages=config["electricity"]["voltages"],
            transformers=config["transformers"],
            snapshots=config["snapshots"],
            links=config["links"],
            lines=config["lines"],
            hvdc_as_lines=config["electricity"]["hvdc_as_lines"],
            countries=config["countries"],
            base_network=config["base_network"],
        input:
            osm_buses="resources/" + RDIR + "base_network/all_buses_build_network.csv",
            osm_lines="resources/" + RDIR + "base_network/all_lines_build_network.csv",
            osm_converters="resources/"
            + RDIR
            + "base_network/all_converters_build_network.csv",
            osm_transformers="resources/"
            + RDIR
            + "base_network/all_transformers_build_network.csv",
            country_shapes="resources/shapes/microgrid_shapes.geojson",
            offshore_shapes=pypsaearth("resources/shapes/offshore_shapes.geojson"),
        output:
            "networks/" + RDIR + "base.nc",
        log:
            "logs/" + RDIR + "base_network.log",
        benchmark:
            "benchmarks/" + RDIR + "base_network"
        threads: 1
        resources:
            mem_mb=500,
        script:
            pypsaearth("scripts/base_network.py")

    rule dist_build_bus_regions:
        params:
            alternative_clustering=config["cluster_options"]["alternative_clustering"],
            crs=config["crs"],
            countries=config["countries"],
        input:
            country_shapes="resources/shapes/microgrid_shapes.geojson",
            offshore_shapes=pypsaearth("resources/shapes/offshore_shapes.geojson"),
            base_network="networks/" + RDIR + "base.nc",
            #gadm_shapes="resources/" + RDIR + "shapes/MAR2.geojson",
            #using this line instead of the following will test updated gadm shapes for MA.
            #To use: downlaod file from the google drive and place it in resources/" + RDIR + "shapes/
            #Link: https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM
            gadm_shapes=pypsaearth("resources/" + RDIR + "shapes/gadm_shapes.geojson"),
        output:
            regions_onshore="resources/" + RDIR + "bus_regions/regions_onshore.geojson",
            regions_offshore="resources/"
            + RDIR
            + "bus_regions/regions_offshore.geojson",
        log:
            "logs/" + RDIR + "build_bus_regions.log",
        benchmark:
            "benchmarks/" + RDIR + "build_bus_regions"
        threads: 1
        resources:
            mem_mb=1000,
        script:
            pypsaearth("scripts/build_bus_regions.py")

    rule dist_filter_data:
        input:
            **{
                f"profile_{tech}": f"resources/renewable_profiles/profile_{tech}.nc"
                for tech in config["tech_modelling"]["general_vre"]
            },
            base_network="networks/base.nc",
            raw_lines="resources/osm/clean/all_clean_lines.geojson",
            shape="resources/shapes/microgrid_shapes.geojson",
        output:
            base_update="networks/" + RDIR + "base_update.nc",
        log:
            "logs/" + RDIR + "dist_filter_data.log",
        benchmark:
            "benchmarks/" + RDIR + "dist_filter_data"
        threads: 1
        resources:
            mem_mb=500,
        script:
            "../scripts/dist_filter_data.py"


rule dist_build_renewable_profiles:
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
        gebco=pypsaearth("data/gebco/GEBCO_2025_sub_ice.nc"),
        country_shapes="resources/shapes/microgrid_shapes.geojson",
        offshore_shapes=pypsaearth("resources/shapes/offshore_shapes.geojson"),
        hydro_capacities="pypsa-earth/data/hydro_capacities.csv",
        eia_hydro_generation="pypsa-earth/data/eia_hydro_annual_generation.csv",
        powerplants="resources/powerplants.csv",
        regions=(
            (
                lambda w: (
                    ("resources/" + RDIR + "bus_regions/regions_onshore.geojson")
                    if w.technology in ("onwind", "solar", "hydro", "csp")
                    else ("resources/" + RDIR + "bus_regions/regions_offshore.geojson")
                )
            )
            if config.get("mode") == "brown_field"
            else "resources/shapes/microgrid_bus_shapes.geojson"
        ),
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


rule dist_add_electricity:
    params:
        mode=config["mode"],
    input:
        **{
            f"profile_{tech}": f"resources/renewable_profiles/profile_{tech}.nc"
            for tech in config["tech_modelling"]["general_vre"]
        },
        create_network=(
            "networks/base_update.nc"
            if config.get("mode") == "brown_field"
            else "networks/base.nc"
        ),
        tech_costs=COSTS,
        load_file="resources/demand/microgrid_load.csv",
        powerplants="resources/powerplants.csv",
    output:
        "networks/elec.nc",
    log:
        "logs/dist_add_electricity.log",
    benchmark:
        "benchmarks/dist_add_electricity"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "../scripts/dist_add_electricity.py"


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


rule dist_solve_network:
    input:
        "networks/elec.nc",
    output:
        "networks/results/elec.nc",
    log:
        "logs/dist_solve_network.log",
    benchmark:
        "benchmarks/dist_solve_network"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "../scripts/dist_solve_network.py"


# if config["enable"].get("retrieve_databundle", True):

#    bundles_to_download = get_best_bundles_in_snakemake(config)

#    rule retrieve_databundle_light:
#        params:
#            bundles_to_download=bundles_to_download,
#            hydrobasins_level=config["renewable"]["hydro"]["hydrobasins_level"],
#        output:  #expand(directory('{file}') if isdir('{file}') else '{file}', file=datafiles)
#            expand(
#                "{file}", file=datafiles_retrivedatabundle(config, bundles_to_download)
#            ),
#            directory("data/landcover"),
#        log:
#            "logs/" + RDIR + "retrieve_databundle.log",
#        benchmark:
#            "benchmarks/" + RDIR + "retrieve_databundle_light"
#        script:
#            "../pypsa-earth/scripts/retrieve_databundle_light.py"

# rule build_shapes:
#    params:
#        build_shape_options=config["build_shape_options"],
#        crs=config["crs"],
#        countries=config["countries"],
#        subregion=config["subregion"],
#    input:
#        # naturalearth='data/bundle/naturalearth/ne_10m_admin_0_countries.shp',
#        # eez='data/bundle/eez/World_EEZ_v8_2014.shp',
#        # nuts3='data/bundle/NUTS_2013_60M_SH/data/NUTS_RG_60M_2013.shp',
#        # nuts3pop='data/bundle/nama_10r_3popgdp.tsv.gz',
#        # nuts3gdp='data/bundle/nama_10r_3gdp.tsv.gz',
#        eez="data/eez/eez_v11.gpkg",
#    output:
#        country_shapes="resources/" + RDIR + "shapes/country_shapes.geojson",
#        offshore_shapes="resources/" + RDIR + "shapes/offshore_shapes.geojson",
#        africa_shape="resources/" + RDIR + "shapes/africa_shape.geojson",
#        gadm_shapes="resources/" + RDIR + "shapes/gadm_shapes.geojson",
#        subregion_shapes="resources/" + RDIR + "shapes/subregion_shapes.geojson",
#    log:
#        "logs/" + RDIR + "build_shapes.log",
#    benchmark:
#        "benchmarks/" + RDIR + "build_shapes"
#    threads: 1
#    resources:
#        mem_mb=3096,
#    script:
#        "scripts/build_shapes.py"

# if config["enable"].get("build_natura_raster", False):
#    rule build_natura_raster:
#        params:
#            area_crs=config["crs"]["area_crs"],
#            natura=config["natura"],
#            disable_progress=not config["enable"]["progress_bar"],
#        input:
#            shapefiles_land="data/landcover",
#            cutouts=expand(
#                "cutouts/" + "{cutout}.nc",
#                cutout=[c["cutout"] for _, c in config["renewable"].items()],
#            ),
#            country_shapes="resources/" + RDIR + "shapes/country_shapes.geojson",
#            offshore_shapes="resources/" + RDIR + "shapes/offshore_shapes.geojson",
#        output:
#            "resources/" + RDIR + "natura.tiff",
#        log:
#            "logs/" + RDIR + "build_natura_raster.log",
#        benchmark:
#            "benchmarks/" + RDIR + "build_natura_raster"
#        script:
#            "scripts/build_natura_raster.py"

# -*- coding: utf-8 -*-
"""
Adds electrical generators, load and storage units to a base network.
Relevant Settings
-----------------
.. code:: yaml
    costs:
        year:
        USD2013_to_EUR2013:
        dicountrate:
    electricity:
        max_hours:
        conventional_carriers:
        extendable_carriers:
    tech_modelling:
        general_vre:
        storage_techs:
        load_carries:

Inputs
------
- ``data/costs.csv``: The database of cost assumptions for all included technologies for specific years from various sources; e.g. discount rate, lifetime, investment (CAPEX), fixed operation and maintenance (FOM), variable operation and maintenance (VOM), fuel costs, efficiency, carbon-dioxide intensity.
- ``resources/powerplants.csv``: confer :ref:`powerplants`
- ``resources/profile_{}.nc``: all technologies in ``config["renewables"].keys()``, confer :ref:`renewableprofiles`
- ``resources/demand/microgrid_load.csv``: microgrid electric demand 
- ``networks/base.nc``: confer :ref:`base`
     
Outputs
-------
- ``networks/elec.nc``:

Description
-----------
The rule :mod:`add_electricity` takes as input the network generated in the rule "create_network" and adds to it both renewable and conventional generation, storage units and load, resulting in a network that is stored in ``networks/elec.nc``. 

"""


import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import powerplantmatching as pm
import pypsa
import xarray as xr
from _helpers_dist import configure_logging, sets_path_to_root

idx = pd.IndexSlice


def calculate_annuity(n, r):
    """
    Calculate the annuity factor for an asset with lifetime n years and
    discount rate of r, e.g. annuity(20, 0.05) * 20 = 1.6
    """
    if isinstance(r, pd.Series):
        return pd.Series(1 / n, index=r.index).where(
            r == 0, r / (1.0 - 1.0 / (1.0 + r) ** n)
        )
    elif r > 0:
        return r / (1.0 - 1.0 / (1.0 + r) ** n)
    else:
        return 1 / n


def _add_missing_carriers_from_costs(n, costs, carriers):
    missing_carriers = pd.Index(carriers).difference(n.carriers.index)
    if missing_carriers.empty:
        return

    emissions_cols = (
        costs.columns.to_series().loc[lambda s: s.str.endswith("_emissions")].values
    )
    suptechs = missing_carriers.str.split("-").str[0]
    emissions = costs.loc[suptechs, emissions_cols].fillna(0.0)
    emissions.index = missing_carriers
    n.import_components_from_dataframe(emissions, "Carrier")


# Last line of costs.csv file is totally invented, it should be reviewed.


def load_costs(tech_costs, config, elec_config, Nyears=1):
    """
    set all asset costs and other parameters
    """
    costs = pd.read_csv(tech_costs, index_col=list(range(3))).sort_index()

    # correct units to MW and EUR
    costs.loc[costs.unit.str.contains("/kW"), "value"] *= 1e3
    costs.loc[costs.unit.str.contains("USD"), "value"] *= config["USD2013_to_EUR2013"]

    costs = (
        costs.loc[idx[:, config["year"], :], "value"]
        .unstack(level=2)
        .groupby("technology")
        .sum(min_count=1)
    )

    costs = costs.fillna(
        {
            "CO2 intensity": 0,
            "FOM": 0,
            "VOM": 0,
            "discount rate": config["discountrate"],
            "efficiency": 1,
            "fuel": 0,
            "investment": 0,
            "lifetime": 25,
        }
    )

    costs["capital_cost"] = (
        (
            calculate_annuity(costs["lifetime"], costs["discount rate"])
            + costs["FOM"] / 100.0
        )
        * costs["investment"]
        * Nyears
    )

    costs.at["OCGT", "fuel"] = costs.at["gas", "fuel"]
    costs.at["CCGT", "fuel"] = costs.at["gas", "fuel"]

    costs["marginal_cost"] = costs["VOM"] + costs["fuel"] / costs["efficiency"]

    costs = costs.rename(columns={"CO2 intensity": "co2_emissions"})

    costs.at["OCGT", "co2_emissions"] = costs.at["gas", "co2_emissions"]
    costs.at["CCGT", "co2_emissions"] = costs.at["gas", "co2_emissions"]

    costs.at["solar", "capital_cost"] = 0.5 * (
        costs.at["solar-rooftop", "capital_cost"]
        + costs.at["solar-utility", "capital_cost"]
    )

    def costs_for_storage(store, link1, link2=None, max_hours=1.0):
        capital_cost = link1["capital_cost"] + max_hours * store["capital_cost"]
        if link2 is not None:
            capital_cost += link2["capital_cost"]
        return pd.Series(
            dict(capital_cost=capital_cost, marginal_cost=0.0, co2_emissions=0.0)
        )

    max_hours = elec_config["max_hours"]
    costs.loc["battery"] = costs_for_storage(
        costs.loc[
            "lithium"
        ],  # line 119 in file costs.csv' which was battery storage was modified into lithium (same values left)
        costs.loc["battery inverter"],
        max_hours=max_hours["battery"],
    )
    max_hours = elec_config["max_hours"]
    costs.loc["battery"] = costs_for_storage(
        costs.loc[
            "lead acid"
        ],  # line 120 in file 'costs.csv' which was battery storage was modified into lithium (same values left)
        costs.loc["battery inverter"],
        max_hours=max_hours["battery"],
    )

    costs.loc["H2"] = costs_for_storage(
        costs.loc["hydrogen storage"],
        costs.loc["fuel cell"],
        costs.loc["electrolysis"],
        max_hours=max_hours["H2"],
    )

    for attr in ("marginal_cost", "capital_cost"):
        overwrites = config.get(attr)
        if overwrites is not None:
            overwrites = pd.Series(overwrites)
            costs.loc[overwrites.index, attr] = overwrites

    return costs


def add_bus_at_center(n, number_microgrids):
    """
    Adds a new bus to each network at the center of the existing buses.
    """

    #Identify the microgrids 
    number_microgrids = len(number_microgrids.keys())
    microgrid_ids = [f"microgrid_{i+1}" for i in range(number_microgrids)]

    # Iterate over each microgrid
    for microgrid_id in microgrid_ids:
        # Select the buses belonging to this microgrid
        microgrid_buses = n.buses.loc[
            n.buses.index.str.contains(f"^{microgrid_id}_bus_")
        ]

        # Compute the centroid of the microgrid buses
        center_x = np.mean(microgrid_buses["x"].values)
        center_y = np.mean(microgrid_buses["y"].values)

        # Create a new bus at the centroid
        center_bus_name = f"new_bus_{microgrid_id}"
        n.add("Bus", center_bus_name, x=center_x, y=center_y, v_nom=0.220)


def attach_wind_and_solar(
    n, costs, input_profiles, tech_modelling, extendable_carriers
):
    """
    This function adds wind and solar generators with the time series "profile_{tech}" to the power network

    """

    # Add any missing carriers from the costs data to the tech_modelling variable
    _add_missing_carriers_from_costs(n, costs, tech_modelling)

    # Get the index of the buses in the power network
    buses_i = n.buses.index

    for tech in tech_modelling:
        # Open the dataset for the current technology from the input_profiles
        with xr.open_dataset(getattr(snakemake.input, "profile_" + tech)) as ds:
            # If the dataset's "bus" index is empty, skip to the next technology
            if ds.indexes["bus"].empty:
                continue

            suptech = tech.split("-", 2)[0]

            # Add the wind and solar generators to the power network
            n.madd(
                "Generator",
                ds.indexes["bus"],
                " " + tech,
                bus=buses_i,
                carrier=tech,
                p_nom_extendable=tech in extendable_carriers["Generator"],
                p_nom_max=ds["p_nom_max"].to_pandas(),  # look at the config
                weight=ds["weight"].to_pandas(),
                marginal_cost=costs.at[suptech, "marginal_cost"],
                capital_cost=costs.at[tech, "capital_cost"],
                efficiency=costs.at[suptech, "efficiency"],
                p_set=ds["profile"]
                .transpose("time", "bus")
                .to_pandas()
                .reindex(n.snapshots),
                p_max_pu=ds["profile"]
                .transpose("time", "bus")
                .to_pandas()
                .reindex(n.snapshots),
            )


def load_powerplants(ppl_fn):
    carrier_dict = {
        "ocgt": "OCGT",
        "ccgt": "CCGT",
        "bioenergy": "biomass",
        "ccgt, thermal": "CCGT",
        "hard coal": "coal",
        # "oil" : "diesel" #This is something that could be done
    }

    return (
        pd.read_csv(ppl_fn, index_col=0, dtype={"bus": "str"})
        .powerplant.to_pypsa_names()
        .powerplant.convert_country_to_alpha2()
        .rename(columns=str.lower)
        .drop(columns=["efficiency"])
        .replace({"carrier": carrier_dict})
    )


def attach_conventional_generators(
    n,
    costs,
    ppl,
    conventional_carriers,
    extendable_carriers,
    conventional_config,
    conventional_inputs,
):
    # Create a set of all conventional and extendable carriers
    carriers = set(conventional_carriers) | set(extendable_carriers["Generator"])

    # Add any missing carriers from the costs data to the "carriers" variable
    _add_missing_carriers_from_costs(n, costs, carriers)

    # Filter the ppl dataframe to only include the relevant carriers
    ppl = (
        ppl.query("carrier in @carriers")
        .join(costs, on="carrier", rsuffix="_r")
        .rename(index=lambda s: "C" + str(s))
    )
    ppl["efficiency"] = ppl.efficiency.fillna(ppl.efficiency)

    # Get the index of the buses in the power network
    buses_i = n.buses.index

    # Add conventional generators to each bus in the power network (one for microgrid)

    n.madd(
        "Generator",
        ppl.index,
        carrier=ppl.carrier,
        bus=buses_i,
        p_nom_min=ppl.p_nom.where(ppl.carrier.isin(conventional_carriers), 0),
        p_nom=ppl.p_nom.where(ppl.carrier.isin(conventional_carriers), 0),
        p_nom_extendable=ppl.carrier.isin(extendable_carriers["Generator"]),
        efficiency=ppl.efficiency,
        marginal_cost=ppl.marginal_cost,
        capital_cost=ppl.capital_cost,
        build_year=ppl.datein.fillna(0).astype(int),
        lifetime=(ppl.dateout - ppl.datein).fillna(np.inf),
    )

    for carrier in conventional_config:
        # Generators with technology affected
        idx = n.generators.query("carrier == @carrier").index

        for attr in list(set(conventional_config[carrier]) & set(n.generators)):
            values = conventional_config[carrier][attr]

            if f"conventional_{carrier}_{attr}" in conventional_inputs:
                # Values affecting generators of technology k country-specific
                # First map generator buses to countries; then map countries to p_max_pu
                values = pd.read_csv(values, index_col=0).iloc[:, 0]
                bus_values = n.buses.country.map(values)
                n.generators[attr].update(
                    n.generators.loc[idx].bus.map(bus_values).dropna()
                )
            else:
                # Single value affecting all generators of technology k indiscriminately of country
                n.generators.loc[idx, attr] = values


def attach_storageunits(n, costs, technologies, extendable_carriers):
    """
    This function adds different technologies of storage units to the power network

    """

    elec_opts = snakemake.config["electricity"]
    max_hours = elec_opts["max_hours"]

    lookup_store = {"H2": "electrolysis", "battery": "battery inverter"}
    lookup_dispatch = {"H2": "fuel cell", "battery": "battery inverter"}

    buses_i = n.buses.index

    # Iterate through each storage technology
    for tech in technologies:
        # Add the storage units to the power network
        n.madd(
            "StorageUnit",
            buses_i,
            " " + tech,
            bus=buses_i,
            carrier=tech,
            p_nom_extendable=True,
            capital_cost=costs.at[tech, "capital_cost"],
            marginal_cost=costs.at[tech, "marginal_cost"],
            efficiency_store=costs.at[
                lookup_store["battery"], "efficiency"
            ],  # Lead_acid and lithium have the same value
            efficiency_dispatch=costs.at[
                lookup_dispatch["battery"], "efficiency"
            ],  # Lead_acid and lithium have the same value
            max_hours=max_hours["battery"],  # Lead_acid and lithium have the same value
            cyclic_state_of_charge=True,
        )


def attach_load(n, load_file, microgrids_list, tech_modelling):
    # Upload the load csv file
    load = pd.read_csv(load_file).set_index([n.snapshots])

    # Create an index for the loads
    index = load.columns

    # Get the index of the buses in the power network
    buses_i = n.buses.index

    # Add the load to the power network
    n.madd("Load", index, bus=buses_i, carrier="AC", p_set=load)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_dist import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("add_electricity")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.create_network)
    Nyears = n.snapshot_weightings.objective.sum() / 8760.0

    load_file = snakemake.input["load_file"]

    ppl = load_powerplants(snakemake.input.powerplants)

    costs = load_costs(
        snakemake.input.tech_costs,
        snakemake.config["costs"],
        snakemake.config["electricity"],
        Nyears,
    )

    add_bus_at_center(n, snakemake.config["microgrids_list"])

    attach_wind_and_solar(
        n,
        costs,
        snakemake.input,
        snakemake.config["tech_modelling"]["general_vre"],
        snakemake.config["electricity"]["extendable_carriers"],
    )

    conventional_inputs = {
        k: v for k, v in snakemake.input.items() if k.startswith("conventional_")
    }

    attach_conventional_generators(
        n,
        costs,
        ppl,
        snakemake.config["electricity"]["conventional_carriers"],
        snakemake.config["electricity"]["extendable_carriers"],
        snakemake.config.get("conventional", {}),
        conventional_inputs,
    )

    attach_storageunits(
        n,
        costs,
        snakemake.config["tech_modelling"]["storage_techs"],
        snakemake.config["electricity"]["extendable_carriers"],
    )

    attach_load(
        n,
        load_file,
        snakemake.config["microgrids_list"],
        snakemake.config["tech_modelling"]["load_carriers"],
    )

    n.export_to_netcdf(snakemake.output[0])

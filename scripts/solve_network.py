# Optimization
import os
from _helpers import configure_logging, sets_path_to_root
import numpy as np
import pypsa
from pypsa.linopf import ilopf, network_lopf

"""
Solves linear optimal power flow for a network iteratively.
-----------------
.. code:: yaml
    solving:
        tmpdir:
        options:
            formulation:
            clip_p_max_pu:
            load_shedding:
            noisy_costs:
            nhours:
            min_iterations:
            max_iterations:
            skip_iterations:
            track_iterations:
        
Inputs
------
- ``networks/elec.nc
Outputs
-------
- ``networks/results/networks/elec.nc : Solved PyPSA network including optimisation results
  
Description
-----------
Total annual system costs are minimised with PyPSA. The full formulation of the
linear optimal power flow (plus investment planning
is provided in the
`documentation of PyPSA <https://pypsa.readthedocs.io/en/latest/optimal_power_flow.html#linear-optimal-power-flow>`_.
The optimization is based on the ``pyomo=False`` setting in the :func:`network.lopf` function.
Additionally, some extra constraints specified in :mod:`prepare_network` are added.
"""

def prepare_network(n, solve_opts):

    if "clip_p_max_pu" in solve_opts:
        for df in (n.generators_t.p_max_pu, n.storage_units_t.inflow):
            df.where(df > solve_opts["clip_p_max_pu"], other=0.0, inplace=True)

    load_shedding = solve_opts.get("load_shedding")
    if load_shedding:
        n.add("Carrier", "Load")
        buses_i = n.buses.query("carrier == 'AC'").index
        if not np.isscalar(load_shedding):
            load_shedding = 8e3  # Eur/kWh
        # intersect between macroeconomic and surveybased
        # willingness to pay
        # http://journal.frontiersin.org/article/10.3389/fenrg.2015.00055/full)
        # 1e2 is practical relevant, 8e3 good for debugging
        n.madd(
            "Generator",
            buses_i,
            " load",
            bus=buses_i,
            carrier="load",
            sign=1e-3,  # Adjust sign to measure p and p_nom in kW instead of MW
            marginal_cost=load_shedding,
            p_nom=1e9,  # kW
        )

    # if solve_opts.get("noisy_costs"):
    #     for t in n.iterate_components(n.one_port_components):
    #         # TODO: uncomment out to and test noisy_cost (makes solution unique)
    #         # if 'capital_cost' in t.df:
    #         #    t.df['capital_cost'] += 1e1 + 2.*(np.random.random(len(t.df)) - 0.5)
    #         if "marginal_cost" in t.df:
    #             t.df["marginal_cost"] += 1e-2 + 2e-3 * (
    #                 np.random.random(len(t.df)) - 0.5
    #             )

        # for t in n.iterate_components(["Line", "Link"]):
        #     t.df["capital_cost"] += (
        #         1e-1 + 2e-2 * (np.random.random(len(t.df)) - 0.5)
        #     ) * t.df["length"]

    # if solve_opts.get("nhours"):
    #     nhours = solve_opts["nhours"]
    #     n.set_snapshots(n.snapshots[:nhours])
    #     n.snapshot_weightings[:] = 8760.0 / nhours

    return n



def solve_network(n, solver_name):

    network_lopf(
            n,
            solver_name=solver_name
        )

    return n
    
if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake(
            "solve_network"
        )

    configure_logging(snakemake)

    solve_opts = snakemake.config["solving"]["options"]

    n = pypsa.Network(snakemake.input[0])
    n = prepare_network(n, solve_opts)
    n = solve_network(
            n,
            "gurobi",
        )

    n.export_to_netcdf(snakemake.output[0]) 
 
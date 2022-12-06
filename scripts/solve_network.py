# Optimization
import os
from _helpers import configure_logging, sets_path_to_root
import numpy as np
import pypsa
from pypsa.linopf import ilopf, network_lopf

def solve_network(n, solver_name):
    # solver_options = config["solving"]["solver"].copy()
    # solver_name = solver_options.pop("name")

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

    # tmpdir = snakemake.config["solving"].get("tmpdir")
    # if tmpdir is not None:
    #     Path(tmpdir).mkdir(parents=True, exist_ok=True)

    n = pypsa.Network(snakemake.input[0])

    n = solve_network(
            n,
            "gurobi",
        )

    n.export_to_netcdf(snakemake.output[0]) 
 
#Creation of the network 

import pypsa
import pandas as pd
import os
from _helpers import configure_logging, sets_path_to_root

def create_network():
    
    n = pypsa.Network()
    n.name = "PyPSA-Distribution"

    n.set_snapshots(pd.date_range(freq="h", **snakemake.config["snapshots"]))
    n.snapshot_weightings[:] *= 8760.0 / n.snapshot_weightings.sum()

    return n

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("create_network")
        
    configure_logging(snakemake)

    n = create_network()
    n.export_to_netcdf(snakemake.output[0])










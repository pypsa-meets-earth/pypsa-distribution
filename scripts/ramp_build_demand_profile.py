# -*- coding: utf-8 -*-
import logging
import os
from pathlib import Path

import pandas as pd
import yaml
from _helpers_dist import configure_logging, create_logger, read_osm_config
from ramp import Appliance, UseCase, User, get_day_type


def create_demand_profile(xlsx_input_path, excel_output_path):
    use_case = UseCase()
    use_case.load(xlsx_input_path)
    # Versione che fa solo una simulazione:
    n_days = 365
    date_start = "2013-01-01"
    use_case.date_start = date_start
    use_case.initialize(num_days=n_days, force=True)
    data = use_case.generate_daily_load_profiles(flat=True)
    data = data / 100
    profile = pd.DataFrame(
        data=data,
        index=pd.date_range(start=date_start, periods=1440 * n_days, freq="T"),
    )

    profile.to_excel(excel_output_path)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_dist import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("ramp_build_demand_profile", user_type="School")
        sets_path_to_root("pypsa-distribution")
    configure_logging(snakemake)

    create_demand_profile(
        snakemake.input["user_description"], snakemake.output["profile_results"]
    )

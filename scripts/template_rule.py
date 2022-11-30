# -*- coding: utf-8 -*-
import os

from _helpers import configure_logging, sets_path_to_root


def main_template_rule(options):
    return 0


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("template_rule")

    configure_logging(snakemake)

    # load default crs
    config = snakemake.config
    input = snakemake.input
    output = snakemake.output

    main_template_rule(config)

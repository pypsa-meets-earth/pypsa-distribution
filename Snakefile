import sys

sys.path.append("./scripts")
sys.path.append("./pypsa-earth/scripts")

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.io import expand
from build_test_configs import create_test_config
from _helpers import create_country_list
from os.path import normpath, exists, isdir
from shutil import copyfile
from pathlib import Path
from _helpers_dist import create_country_list
import os

import sys

sys.path.append("./scripts")

from os.path import normpath, exists, isdir

HTTP = HTTPRemoteProvider()


COSTS = "data/costs.csv"
PROFILE = "data/sample_profile.csv"

PYPSAEARTH_FOLDER = "pypsa-earth"

if not config.get("disable_subworkflow", False):

    subworkflow pypsaearth:
        workdir:
            PYPSAEARTH_FOLDER
        snakefile:
            PYPSAEARTH_FOLDER + "/Snakefile"
        configfile:
            "configs/config.pypsa-earth.yaml"


if config.get("disable_subworkflow", False):

    def pypsaearth(path):
        return PYPSAEARTH_FOLDER + "/" + path


configfile: "configs/config.pypsa-earth.yaml"


configfile: "configs/config.distribution.yaml"


if exists("config.yaml"):

    configfile: "config.yaml"


config["countries"] = create_country_list(config["countries"])

run = config.get("run", {})
RDIR = run["name"] + "/" if run.get("name") else ""
countries = config["countries"]

ATLITE_NPROCESSES = config["atlite"].get("nprocesses", 5)


wildcard_constraints:
    ll="[a-z0-9\.]+",
    simpl="[a-zA-Z0-9]*|all",
    clusters="[0-9]+m?|all",
    opts="[-+a-zA-Z0-9]*",
    sopts="[-+a-zA-Z0-9\.\s]*",
    user_type="[a-zA-Z0-9]*",


if not config.get("disable_subworkflow", False):

    subworkflow pypsaearth:
        workdir:
            PYPSAEARTH_FOLDER
        snakefile:
            PYPSAEARTH_FOLDER + "/Snakefile"
        configfile:
            "./configs/config.pypsa-earth.yaml"


if config.get("disable_subworkflow", False):

    def pypsaearth(path):
        return PYPSAEARTH_FOLDER + "/" + path


if not config["enable"].get("disable_distribution_workflow"):

    include: "rules/pypsa_distribution.smk"


# rule clean:
#     run:
#         shell("snakemake -j 1 solve_network --delete-all-output")

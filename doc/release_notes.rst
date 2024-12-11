..
  SPDX-FileCopyrightText: 2021 The PyPSA-Earth Authors

  SPDX-License-Identifier: CC-BY-4.0

##########################################
Release Notes
##########################################


Upcoming Release
================

* The energy demand time series can now be determined using a new methodology based on social inputs that integrates with the RAMP tool. `PR #55 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/55>`__

* Automated downloading of buildings within the microgrid is now supported through the new download_osm_data rule. `PR #52 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/52>`__ and `PR #56 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/56>`__



Release Process
===============

* Checkout a new release branch ``git checkout -b release-v0.x.x``.

* Finalise release notes at ``doc/release_notes.rst``.

* Update ``envs/environment.fixed.yaml`` via
  ``conda env export -n pypsa-earth -f envs/environment.fixed.yaml --no-builds``
  from an up-to-date `pypsa-earth` environment.

* Update version number in ``doc/conf.py`` and ``*config.*.yaml``.

* Open, review and merge pull request for branch ``release-v0.x.x``.
  Make sure to close issues and PRs or the release milestone with it (e.g. closes #X).

* Tag a release on Github via ``git tag v0.x.x``, ``git push``, ``git push --tags``. Include release notes in the tag message.

* Upload code to `zenodo code repository <https://doi.org>`_ with `GPLv3 license <https://www.gnu.org/licenses/gpl-3.0.en.html>`_.

* Send announcement on the `PyPSA-Earth Discord channel <https://discord.gg/AnuJBk23FU>`_.

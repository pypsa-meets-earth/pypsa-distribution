..
  SPDX-FileCopyrightText: 2021 The PyPSA-Earth Authors

  SPDX-License-Identifier: CC-BY-4.0

##########################################
Release Notes
##########################################


Upcoming Release
================

**New Features and Major Changes**

* 

**Minor Changes and bug-fixing**

* 

Version 0.0.2
=================

**New Features and Major Changes**

* Procedure to create brownfield network using OSM data. `PR #75 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/75>`__

* Download network data for distribution systems using OSM. `PR #73 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/73>`__

* Enable interconnection of multiple microgrids with kmeans. `PR #69 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/69>`__

* Enable interconnection of multiple microgrids. `PR #70 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/70>`__

* Implementation of documentation `PR #61 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/61>`__ and `PR #62 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/62>`__

* Introduce generation bus `PR #59 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/59>`__

**Minor Changes and bug-fixing**

* Improve README `PR #65 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/65>`__

* Bug-fixing build demand `PR #63 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/63>`__

* Prepare v0.0.2. `PR #77 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/77>`__


Version 0.0.1
=================

**New Features and Major Changes**

* The energy demand time series can now be determined using a new methodology based on social inputs that integrates with the RAMP tool. `PR #55 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/55>`__

* Automated downloading of buildings within the microgrid is now supported through the new download_osm_data rule. `PR #52 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/52>`__ and `PR #56 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/56>`__

* Introduce demand estimation by ramp `PR #45 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/45>`__

* First draft to download OSM buildings using eart-osm `PR #41 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/41>`__

* Introduce gadm shape `PR #4 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/4>`__

* Automatize build renewable production `PR #2 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/2>`__

* Workflow execution `PR #1 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/1>`__

**Minor Changes and bug-fixing**

* Bug-fix multi-microgrid `PR #58 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/58>`__

* Add rampdemand installation guide `PR #49 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/49>`__

* Add README and PR template `PR #38 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/38>`__ and `PR #37 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/37>`__

* Add option to disable subworkflow `PR #36 <https://github.com/pypsa-meets-earth/pypsa-distribution/pull/36>`__


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

# MONAN-v0.1.0

## Model for Ocean-laNd-Atmosphere predictioN

MONAN is a community model of the Unified Earth System that has as its main objective, as its name already refers, to be a numerical model that covers all scales, geographical and temporal, of the entire Earth system and its implications. It is "community" because it aggregates efforts from several brazilian national institutions such as universities, research centers, operational centers and various authorities in the area of meteorology, environment, oceans and others. It also can receive support from international centers and universities as well as support from the private sector.

The MONAN Model is managed by a scientific committee appointed by INPE's director and has its initial version structure (0.1.0) based on the dynamic core of the MPAS 8.0.1 Model. part of the physics used by MONAN is obtained from the MPAS model and another part obtained from other sources or developed by the community. The MPAS model can be found at the link [GitHub - MPAS-Dev/MPAS-Model: Repository for MPAS models and shared framework releases.](https://github.com/MPAS-Dev/MPAS-Model)

MPAS-v8.0.1
====

The Model for Prediction Across Scales (MPAS) is a collaborative project for
developing atmosphere, ocean, and other earth-system simulation components for
use in climate, regional climate, and weather studies. The primary development
partners are the climate modeling group at Los Alamos National Laboratory
(COSIM) and the National Center for Atmospheric Research. Both primary
partners are responsible for the MPAS framework, operators, and tools common to
the applications; LANL has primary responsibility for the ocean model, and NCAR
has primary responsibility for the atmospheric model.

The MPAS framework facilitates the rapid development and prototyping of models
by providing infrastructure typically required by model developers, including
high-level data types, communication routines, and I/O routines. By using MPAS,
developers can leverage pre-existing code and focus more on development of
their model.

BUILDING
========

This README is provided as a brief introduction to the MPAS framework. It does
not provide details about each specific model, nor does it provide building
instructions.

For information about building and running each core, please refer to each
core's user's guide, which can be found at the following web sites:

[MPAS-Atmosphere](http://mpas-dev.github.io/atmosphere/atmosphere_download.html)

[MPAS-Albany Land Ice](http://mpas-dev.github.io/land_ice/download.html)

[MPAS-Ocean](http://mpas-dev.github.io/ocean/releases.html)

[MPAS-Seaice](http://mpas-dev.github.io/sea_ice/releases.html)

Code Layout
----------

Within the MPAS repository, code is laid out as follows. Sub-directories are
only described below the src directory.

    MPAS-Model
    ├── src
    │   ├── driver -- Main driver for MPAS in stand-alone mode (Shared)
    │   ├── external -- External software for MPAS (Shared)
    │   ├── framework -- MPAS Framework (Includes DDT Descriptions, and shared routines. Shared)
    │   ├── operators -- MPAS Opeartors (Includes Operators for MPAS meshes. Shared)
    │   ├── tools -- Empty directory for include files that Registry generates (Shared)
    │   │   ├── registry -- Code for building Registry.xml parser (Shared)
    │   │   └── input_gen -- Code for generating streams and namelist files (Shared)
    │   └── core_* -- Individual model cores.
    │       └── inc -- Empty directory for include files that Registry generates
    ├── testing_and_setup -- Tools for setting up configurations and test cases (Shared)
    └── default_inputs -- Copies of default stream and namelists files (Shared)

Model cores are typically developed independently. For information about
building and running a particular core, please refer to that core's user's
guide.

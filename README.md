[![Build Status](https://travis-ci.org/pandegroup/openmm.svg?branch=master)](https://travis-ci.org/pandegroup/openmm)
[![Anaconda Cloud Badge](https://anaconda.org/omnia/openmm/badges/downloads.svg)](https://anaconda.org/omnia/openmm)

# Tinker-OpenMM: Version of OpenMM for Use with Tinker

<H2><B><I>NOTE: This OpenMM version for use as part of Tinker-OpenMM is now deprecated. The current Tinker release on the TinkerTools Github site supports Tinker-OpenMM using the current release of OpenMM available from https://github.com/openmm. The modified OpenMM located here should only be used and will only work correctly with Tinker releases prior to Tinker 8.8.</I></B></H2>

<H2><B>Introduction</B></H2>

[OpenMM](http://openmm.org) is a toolkit for molecular simulation. It can be used either as a stand-alone application for running simulations, or as a library you call from your own code. It provides a combination of flexibility (through custom forces and integrators), openness, and high performance on recent GPUs. This Github repository differs from the Main OpenMM release as it contains a softcore VDW for lambda based free energy simulations, and other changes ported from the Tinker molecular modeling software. This version of OpenMM can be linked against Tinker to provide executables (dynamic_omm.x, etc.) for runing simulations on GPUs using Tinker as the front-end interface. To build and install this software, which we call Tinker-OpenMM, see instructions at https://biomol.bme.utexas.edu/tinkergpu/index.php?title=Tinkergpu:Tinker-tut

Need Help? Check out the [documentation](http://docs.openmm.org/) and [discussion forums](https://simtk.org/forums/viewforum.php?f=161).

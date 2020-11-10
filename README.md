[![Build Status](https://travis-ci.org/pandegroup/openmm.svg?branch=master)](https://travis-ci.org/pandegroup/openmm)
[![Anaconda Cloud Badge](https://anaconda.org/omnia/openmm/badges/downloads.svg)](https://anaconda.org/omnia/openmm)

# Tinker-OpenMM: Version of OpenMM for Use with Tinker

<H2><B><I>NOTE: This OpenMM version for use as part of Tinker-OpenMM is now deprecated. The current Tinker release on the TinkerTools Github site supports Tinker-OpenMM using the current release of OpenMM available from https://github.com/openmm. The modified OpenMM located here should only be used and will only work correctly with Tinker release 8.7 and earlier.</I></B></H2>

<H2><B>Introduction</B></H2>

[OpenMM](http://openmm.org) is a toolkit for molecular simulation. It can be used either as a stand-alone application for running simulations, or as a library you call from your own code. This Github repository differs from the earlier canonical OpenMM releases as it contains changes ported from the Tinker molecular modeling software to allow use of some Tinker methods and to implement free energy calculations. This version of OpenMM can be linked against Tinker to provide executables (dynamic_omm.x, analyze_omm.x and bar_omm.x) for running simulations on GPUs using Tinker as the front-end interface. To build and install this software, called Tinker-OpenMM, see instructions at https://biomol.bme.utexas.edu/tinkergpu/index.php?title=Tinkergpu:Tinker-tut

Need Help? Check out the [documentation](http://docs.openmm.org/) and [discussion forums](https://simtk.org/forums/viewforum.php?f=161).

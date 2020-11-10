[![Build Status](https://travis-ci.org/pandegroup/openmm.svg?branch=master)](https://travis-ci.org/pandegroup/openmm)
[![Anaconda Cloud Badge](https://anaconda.org/omnia/openmm/badges/downloads.svg)](https://anaconda.org/omnia/openmm)

# Tinker-OpenMM: Version of OpenMM for Use with Tinker

<H2><B><I>NOTE: This OpenMM version for use as part of Tinker-OpenMM is now deprecated. The current Tinker release on the TinkerTools Github site supports Tinker-OpenMM using the current release of OpenMM available from https://github.com/openmm. The modified OpenMM located here should only be used and will only work correctly with Tinker release 8.7 and earlier.</I></B></H2>

<H2><B>Introduction</B></H2>

Tinker-OpenMM provides a frontend to [OpenMM](http://openmm.org), a toolkit for molecular simulation. As such Tinker-OpenMM allows use of Tinker files, parameters and algorithms to perform modeling calculations and molecular dynamics simulations with OpenMM as the back-end compute engine. This repository differs from the earlier canonical OpenMM releases in that it contains changes ported from the Tinker molecular modeling software to allow use of selected Tinker methods and to implement free energy calculations. The OpenMM version provided here is to be linked against Tinker 8.7 or earlier to provide Tinker-OpenMM executables (dynamic_omm.x, analyze_omm.x and bar_omm.x). To build and install this software, see the instructions and documentation available at https://biomol.bme.utexas.edu/tinkergpu/index.php?title=Tinkergpu:Tinker-tut.

Need help building OpenMM? Check out the OpenMM [documentation](http://docs.openmm.org/) and [discussion forums](https://simtk.org/forums/viewforum.php?f=161).

This branch is the preliminary implementation of AMOEBA+ model with the latest charge flux (CF) implementation.

We are currently migrating our implementation into the latest OpenMM.

Date: 11/22/2019

Progress on migrating to the latest OpenMM
- [x] Van der Waals: W-H combining rules now supported in OpenMM
- [ ] Charge tranfer: code optimization using the *nonbonded* force supported by OpenMM
- [ ] AMOEBAplusNonbond: electrostatics and polarization with charge flux

## Reference
1. Liu, C.; Piquemal, J.-P.; Ren, P., Implementation of Geometry Dependent Charge Flux into AMOEBA+ Potential.  *submitted*
1. Liu, C.; Piquemal, J.-P.; Ren, P., AMOEBA+ Classical Potential for Modeling Molecular Interactions. *J. Chem. Theory Comput.* **2019**, 15 (7), 4122-4139(__Vdw, CT, CP, Pol__)
1. Liu, C.; Qi, R.; Wang, Q.; Piquemal, J.-P.; Ren, P., Capturing Many-body Interactions with Classical Dipole Induction Models. *J. Chem. Theory Comput.*, **2017**, 13 (6), 2751-2761 (__Pol__)
1. Rackers, J. A.; Wang, Q.; Liu, C.; Piquemal, J.-P.; Ren, P.; Ponder, J. W., An optimized charge penetration model for use with the AMOEBA force field. *Phys. Chem. Chem. Phys.* **2016**, 19, 276-291 (__CP__)

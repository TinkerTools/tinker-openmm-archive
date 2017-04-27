This folder contains tests for the dual topology features of openMM
Both examples use vdw-lambda= 0.5 and both ele-lambdas as 0 in order to provide confirmation that dual topology simulation is working
If you are not using a version of openMM that supports dual topology, the total energy will be large due to ligand1 ligand2 vdW overlap
Remember when running dual topology free energy simulations to never have both ligands charged at the same time
This is because the electrostatic force doesn't have the capability of setting ligand1-ligand2 interactions to zero force and energy. 

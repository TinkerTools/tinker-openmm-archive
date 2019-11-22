source ~/.openmm.amoebaplus2019.cflux
rm -rf *.arc *.dyn
#../../dynamic_omm.x 001_water.xyz -key liquid.key 1 1 0.001 2 298 
#../../dynamic_omm.x dimer.xyz -key liquid.key 1 1 0.001 2 298 
#../../dynamic_omm.x elevenmer.xyz -key liquid.key 1 1 0.001 2 298 
#../../dynamic_omm.x 10mer.xyz -key liquid.key 1 1 0.001 2 298 
#../../dynamic_omm.x 6mer.xyz -key liquid.key 1 1 0.001 2 298 
#../../dynamic_omm.x 3mer.xyz -key liquid.key 1 1 0.001 2 298 
#../../dynamic_omm.x dimer.xyz -key liquid.key 1 1 0.001 2 298 
../../dynamic_omm.x liquid.xyz -key liquid.key 1 1 0.001 2 298 
#../../dynamic_omm.x 32mer.xyz -key liquid.key 1 1 0.001 2 298 
#~/bin/tinkerAMOEBA+/dynamic_omm liquid.xyz -key liquid.key 1 1 0.001 2 298 ##correct
#../dynamic_omm.x dimer.xyz -key liquid.key 1 1 0.001 4 298 1 
#~/bin/tinkerOMM2019/dynamic_omm dimer.xyz -key liquid.key 1 1 0.001 4 298 1 

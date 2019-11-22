source ~/.bashrc.intel
#../../../Tinker-CPU/analyze.x 001_water.arc -key tinker.key E  
#../../../Tinker-CPU/testgrad.x 001_water.arc -key tinker.key Y N 0.0001  
~/bin/tinkerAMOEBA+/testgrad elevenmer.xyz -key liquid.key Y Y 0.0001  
../../../Tinker-CPU/testgrad.x elevenmer.xyz -key liquid.key Y Y 0.0001  
#../../../Tinker-CPU/testgrad.x liquid.arc -key liquid.key Y N 0.0001  
#../../../Tinker-CPU/testgrad.x liquid.xyz -key liquid.key Y Y 0.0001  
#../../../Tinker-CPU/testgrad.x first.xyz -key liquid.key Y N 0.0001  
#../../../Tinker-CPU/testgrad.x second.xyz -key liquid.key Y N 0.0001  

#~/Softwares/tinkers/tinker8.2-intel13/AMOEBAplus_CFlux/source-CFlux/analyze.x 001_water.arc E 
#~/Softwares/tinkers/tinker8.2-intel13/AMOEBAplus_CFlux/source-CFlux/analyze.x dimer.arc E 

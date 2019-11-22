#AMOEBA
rm -fr *.dyn
source /home/liuchw/.openmm.amoebaplus2019.cflux 
export CUDA_VISIBLE_DEVICES=0
../dynamic_omm.x liquid.xyz 5 1.0 0.001 4 298.15 1 N >GPU_liq.out 


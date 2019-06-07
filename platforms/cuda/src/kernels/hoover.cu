extern "C" __global__ void scaleVelocity(double* __restrict__ Velscale, mixed4* __restrict__ velm, int numAtoms){
	mixed4 velocity;
	for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
	velocity=velm[index];
 	velocity.x*=(mixed)Velscale[0];
	velocity.y*=(mixed)Velscale[0];
	velocity.z*=(mixed)Velscale[0];
        velm[index]=velocity;
     }
}extern "C" __global__ void scalePos(double* __restrict__ Posscale, real4* __restrict__ posq, real4* __restrict__ posqCorrection,int numAtoms){
	real4 pos;
	for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
	pos=posq[index];
	real4 posCorrection=posqCorrection[index];
        mixed4 ActualPos= make_mixed4(pos.x+(mixed) posCorrection.x,pos.y+(mixed) posCorrection.y,pos.z+(mixed) posCorrection.z,pos.w);
 	ActualPos.x*=Posscale[0];
	ActualPos.y*=Posscale[0];
	ActualPos.z*=Posscale[0];
	posq[index]=make_real4((real) ActualPos.x,(real) ActualPos.y,(real) ActualPos.z,pos.w);
        posqCorrection[index]= make_real4( ActualPos.x-(real) ActualPos.x,ActualPos.y-(real) ActualPos.y,ActualPos.z-(real) ActualPos.z,pos.w);        
     }
}

#include <stdio.h>
extern "C" __global__ void printatom(const long long* __restrict__ force,mixed3* slowaccel, mixed3* fastaccel, real4* __restrict__ posq){
	for( int index=0; index<9; index++){
		    const mixed scale = 0.001/(mixed) 0x100000000;
        	printf("%g   %g   %g  %g  %g  %g %g %g %g \n",posq[index].x,posq[index].y,posq[index].z, slowaccel[index].x,slowaccel[index].y,slowaccel[index].z,fastaccel[index].x,fastaccel[index].y,fastaccel[index].z);
	}
}
extern "C" __global__ void printatom2(const long long* __restrict__ force, real4* __restrict__ posq, int paddednum){
        for( int index=0; index<9; index++){
                    const mixed scale = 0.001/(mixed) 0x100000000;
                printf("%g   %g   %g  %g  %g  %g %g %g %g \n",posq[index].x,posq[index].y,posq[index].z, scale*force[index],scale*force[index+paddednum],scale*force[index+2*paddednum]);
        }
}
/**
 * Perform the first step of Verlet integration.
 */
extern "C" __global__ void saveaccel(int numAtoms, int paddedNumAtoms, const mixed* __restrict__ StepSizes, mixed4* __restrict__ velm,const long long* __restrict__ force,  mixed3* fastaccel) {
		    const mixed stepSize = StepSizes[1]*(mixed)0.5;
    const mixed scale = stepSize/(mixed) 0x100000000;
		         for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
					    fastaccel[index].x=scale*force[index]*velm[index].w;
	    fastaccel[index].y=scale*force[index+paddedNumAtoms]*velm[index].w;
	    fastaccel[index].z=scale*force[index+paddedNumAtoms*2]*velm[index].w;	 
}
}
extern "C" __global__ void saveslowaccel(int numAtoms, int paddedNumAtoms, const mixed* __restrict__ StepSizes, mixed4* __restrict__ velm,const long long* __restrict__ force,  mixed3* slowaccel){
                    const mixed stepSize = StepSizes[0]*(mixed)0.5;
    const mixed scale = stepSize/(mixed) 0x100000000;
                         for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
           	 slowaccel[index].x=scale*force[index]*velm[index].w;
           	 slowaccel[index].y=scale*force[index+paddedNumAtoms]*velm[index].w;
           	 slowaccel[index].z=scale*force[index+paddedNumAtoms*2]*velm[index].w;

	}
}

extern "C" __global__ void integrateRESPAPart1(int numAtoms, int paddedNumAtoms, const mixed* __restrict__ StepSizes, mixed4* __restrict__ velm, const long long* __restrict__ force,  mixed3* slowaccel) {
    const mixed stepSize = StepSizes[0]*(mixed)0.5;
    const mixed scale = stepSize/(mixed) 0x100000000;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            velocity.x += slowaccel[index].x;
            velocity.y += slowaccel[index].y;
            velocity.z += slowaccel[index].z;
            velm[index] = velocity;
//	    if( index<9){
//		printf("%g  %g   %g\n", slowaccel[index].x,slowaccel[index].y,slowaccel[index].z);
//	    }
        }
    }
}
extern "C" __global__ void integrateRESPAPart2(int numAtoms, int paddedNumAtoms, const mixed* __restrict__ StepSizes, mixed4* __restrict__ velm, const long long* __restrict__ force, real4* __restrict__ posq, real4*  __restrict__ posqCorrection, mixed3* fastaccel) {
    const mixed stepSize = StepSizes[1]*(mixed)0.5;
    const mixed scale = stepSize/(mixed) 0x100000000;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
        mixed4 velocity = velm[index];
	mixed3 acceleration;        
	if (velocity.w != 0.0) {
            velocity.x += fastaccel[index].x;
            velocity.y += fastaccel[index].y;
            velocity.z += fastaccel[index].z;
            velm[index] = velocity;
        }
       real4 pos=posq[index];
       real4 posCorrection=posqCorrection[index];
       mixed4 ActualPos= make_mixed4(pos.x+(mixed) posCorrection.x,pos.y+(mixed) posCorrection.y,pos.z+(mixed) posCorrection.z,pos.w);
        ActualPos.x+=stepSize*velm[index].x*2.0;
        ActualPos.y+=stepSize*velm[index].y*2.0;
        ActualPos.x+=stepSize*velm[index].z*2.0;
 	posq[index]=make_real4((real) ActualPos.x,(real) ActualPos.y,(real) ActualPos.z,pos.w);
        posqCorrection[index]= make_real4( ActualPos.x-(real) ActualPos.x,ActualPos.y-(real) ActualPos.y,ActualPos.z-(real) ActualPos.z,pos.w);      
      }
}
extern "C" __global__ void integrateRESPAPart3(int numAtoms, int paddedNumAtoms, const mixed* __restrict__ StepSizes, mixed4* __restrict__ velm, const long long* __restrict__ force, mixed3* fastaccel, mixed3* slowaccel) {
    const mixed stepSize = StepSizes[1]*(mixed)0.5;
    const mixed scale = stepSize/(mixed) 0x100000000;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
        mixed4 velocity = velm[index];
	mixed3 acceleration;        
	if (velocity.w != 0.0) {
	    acceleration.x=scale*force[index]*velocity.w;
	    acceleration.y= scale*force[index+paddedNumAtoms]*velocity.w;
	    acceleration.z=scale*force[index+paddedNumAtoms*2]*velocity.w;
            velocity.x += acceleration.x;
            velocity.y += acceleration.y;
            velocity.z += acceleration.z;
            velm[index] = velocity;
	    fastaccel[index]=acceleration;
        }    
      }
}extern "C" __global__ void integrateRESPAPart4(int numAtoms, int paddedNumAtoms, const mixed* __restrict__ StepSizes, mixed4* __restrict__ velm, const long long* __restrict__ force, mixed3* slowaccel) {
  	  const mixed stepSize = StepSizes[0]*(mixed)0.5;
  	  const mixed scale = stepSize/(mixed) 0x100000000;	
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
	mixed4 velocity=velm[index];
       	mixed3 acceleration;
	acceleration.x=scale*force[index]*velocity.w;
	acceleration.y= scale*force[index+paddedNumAtoms]*velocity.w;
	acceleration.z=scale*force[index+paddedNumAtoms*2]*velocity.w;
	slowaccel[index]=acceleration;
        velocity.x += acceleration.x;
        velocity.y += acceleration.y;
        velocity.z += acceleration.z;
	velm[index]=velocity;
    }
}

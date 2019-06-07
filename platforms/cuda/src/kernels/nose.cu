#include <stdio.h>
extern "C" __global__ void printatom(mixed4* __restrict__ velm, real4* __restrict__ posq){
	printf("%g   %g   %g  %g  %g  %g",velm[0].x,velm[0].y,velm[0].z,posq[0].x,posq[0].y,posq[0].z);
}
extern "C" __global__ void scaleVelocity(double* __restrict__ Velscale, mixed4* __restrict__ velm, int numAtoms){
	mixed4 velocity;
	for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
	velocity=velm[index];
 	velocity.x*=(mixed)Velscale[0];
	velocity.y*=(mixed)Velscale[0];
	velocity.z*=(mixed)Velscale[0];
        velm[index]=velocity;
     }
}
//extern "C" __global__ void integratePart1(double dt, int numAtoms,float poly, float eterm2,int paddedNumAtoms,mixed4* __restrict__ velm, real4* __restrict__ posq, const long long* __restrict__ force, float vbar){
extern "C" __global__ void integratePart1(int numAtoms,int paddedNumAtoms,double* __restrict__ dt,double* __restrict__ PolyEterm2,mixed4* __restrict__ velm, real4* __restrict__ posq, real4* posqCorrection,const long long* __restrict__ force){
    const double scale = dt[0]/(mixed) 0x100000000;
  for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
	mixed4 velocity=velm[index];
	double mass= velm[index].w;
	velocity.x+=scale*force[index]*mass*0.5;
        velocity.y+=scale*force[index+paddedNumAtoms]*mass*0.5;
	velocity.z+=scale*force[index+paddedNumAtoms*2]*mass*0.5;
	velocity.w=mass;
        real4 pos=posq[index];
	real4 posCorrection=posqCorrection[index];
	mixed4 ActualPos= make_mixed4(pos.x+(mixed) posCorrection.x,pos.y+(mixed) posCorrection.y,pos.z+(mixed) posCorrection.z,pos.w);
	ActualPos.x=ActualPos.x*PolyEterm2[1]+velocity.x*PolyEterm2[0];
	ActualPos.y=ActualPos.y*PolyEterm2[1]+velocity.y*PolyEterm2[0];
	ActualPos.z=ActualPos.z*PolyEterm2[1]+velocity.z*PolyEterm2[0];
	pos.w=posq[index].w;
	posq[index]=make_real4((real) ActualPos.x,(real) ActualPos.y,(real) ActualPos.z,pos.w);
	posqCorrection[index]= make_real4( ActualPos.x-(real) ActualPos.x,ActualPos.y-(real) ActualPos.y,ActualPos.z-(real) ActualPos.z,pos.w);
	velm[index]=velocity;
  }
}extern "C" __global__ void integratePart2(int numAtoms, int paddedNumAtoms, double* __restrict__  dt,mixed4* __restrict__ velm,const mixed4* __restrict__ posq, const long long* __restrict__ force){
	    const double scale = dt[0]/(mixed) 0x100000000;
            mixed3 acceleration;
	    mixed4 velocity;
            for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
                acceleration.x= scale*force[index]*velm[index].w;
                acceleration.y=scale*force[index+paddedNumAtoms]*velm[index].w;
                acceleration.z=scale*force[index+2*paddedNumAtoms]*velm[index].w;
		velocity.x= velm[index].x+acceleration.x*0.5;
                velocity.y= velm[index].y+acceleration.y*0.5;
                velocity.z= velm[index].z+acceleration.z*0.5;
		velocity.w=velm[index].w;
		velm[index]=velocity;
            }
}

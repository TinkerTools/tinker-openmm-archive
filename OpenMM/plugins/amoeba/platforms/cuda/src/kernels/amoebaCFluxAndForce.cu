
// loop over bond to accumulate charge fluxes
extern "C" __global__ void computeBondCFluxes(const real4* __restrict__ posq, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, real* __restrict__ cfCharges, const real2* bondCFluxParams, const uint2* __restrict__ atomIndicesBond) {
    for (unsigned int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_BONDS; index += blockDim.x*gridDim.x) 
    {
        uint2 atoms = atomIndicesBond[index];
        unsigned int atom1 = atoms.x;
        unsigned int atom2 = atoms.y;
        real4 pos1 = posq[atom1];
        real4 pos2 = posq[atom2];
        real3 delta = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
    
#if APPLY_PERIODIC
        APPLY_PERIODIC_TO_DELTA(delta)
#endif
        real r = SQRT(delta.x*delta.x + delta.y*delta.y + delta.z*delta.z);
        float2 data = bondCFluxParams[index];
        real dr = r - data.x;
        real deltaQ = data.y*dr;
         
        real cfQ1 = (-deltaQ); 
        real cfQ2 = -cfQ1;
        
        atomicAdd(&cfCharges[atom1], cfQ1);
        __threadfence_block();
        atomicAdd(&cfCharges[atom2], cfQ2);
        __threadfence_block();
    }
}


// loop over angle to accumulate charge fluxes
extern "C" __global__ void computeAngleCFluxes(const real4* __restrict__ posq, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, const real3* angleCFluxAngleParams, const real4* angleCFluxBondParams, real* __restrict__ cfCharges, const uint3* __restrict__ atomIndicesAngle) {
    for (unsigned int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_ANGLES; index += blockDim.x*gridDim.x) {
        uint3 atoms0 = atomIndicesAngle[index];
        unsigned int atom1 = atoms0.x;
        unsigned int atom2 = atoms0.y;
        unsigned int atom3 = atoms0.z;
        real4 pos1 = posq[atom1];
        real4 pos2 = posq[atom2];
        real4 pos3 = posq[atom3];
        real3 v0 = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
        real3 v1 = make_real3(pos2.x-pos3.x, pos2.y-pos3.y, pos2.z-pos3.z);
#if APPLY_PERIODIC
        APPLY_PERIODIC_TO_DELTA(v0)
        APPLY_PERIODIC_TO_DELTA(v1)
#endif
        real r21 = v0.x*v0.x + v0.y*v0.y + v0.z*v0.z;
        real r23 = v1.x*v1.x + v1.y*v1.y + v1.z*v1.z;
        real dot = v0.x*v1.x + v0.y*v1.y + v0.z*v1.z;
        real cosine = min(max(dot*RSQRT(r21*r23), (real) -1), (real) 1);
    
        real theta = ACOS(cosine)*RAD_TO_DEG; 

        r21 = SQRT(r21);
        r23 = SQRT(r23);

        float3 angleParams = angleCFluxAngleParams[index];
        real theta0 =  angleParams.x;
        real jtheta1 = angleParams.y;
        real jtheta2 = angleParams.z;
     
        float4 bondParams = angleCFluxBondParams[index];
        real bond1 = bondParams.x;
        real jbp1  = bondParams.y;
        real bond2 = bondParams.z;
        real jbp2  = bondParams.w;
        
        real dQ1 = jbp1*(r23-bond2) + jtheta1*(theta - theta0);
        real dQ3 = jbp2*(r21-bond1) + jtheta2*(theta - theta0);
        real dQ2 = -(dQ1+dQ3);
        
        atomicAdd(&cfCharges[atom1], dQ1);
        __threadfence_block();
        atomicAdd(&cfCharges[atom2], dQ2);
        __threadfence_block();
        atomicAdd(&cfCharges[atom3], dQ3);
        __threadfence_block();
    }
}


//loop over each bond to calculate chain rule term force 
extern "C" __global__ void computeCFluxBondForce(unsigned long long* __restrict__ forceBuffers, const real4* __restrict__ posq, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, float* __restrict__ cfPotentials, const float2* bondCFluxParams, const uint2* __restrict__ atomIndicesBond) {
    for (unsigned int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_BONDS; index += blockDim.x*gridDim.x) 
    {
        uint2 atoms = atomIndicesBond[index];
        unsigned int atom1 = atoms.x;
        unsigned int atom2 = atoms.y;
        real4 pos1 = posq[atom1];
        real4 pos2 = posq[atom2];
        // b - a
        real3 delta = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
    
#if APPLY_PERIODIC
        APPLY_PERIODIC_TO_DELTA(delta)
#endif
        real r = SQRT(delta.x*delta.x + delta.y*delta.y + delta.z*delta.z);
        float2 data = bondCFluxParams[index];
        real deltaPot = cfPotentials[atom1]-cfPotentials[atom2];
        real3 ddqdr = make_real3(data.y*delta.x/r, data.y*delta.y/r,data.y*delta.z/r);
        real3 force1 = make_real3(-ddqdr.x*deltaPot, -ddqdr.y*deltaPot, -ddqdr.z*deltaPot);
        real3 force2 = -force1;
        atomicAdd(&forceBuffers[atom1],                    static_cast<unsigned long long>((long long) (force1.x*0x100000000)));
        atomicAdd(&forceBuffers[atom1+PADDED_NUM_ATOMS],   static_cast<unsigned long long>((long long) (force1.y*0x100000000)));
        atomicAdd(&forceBuffers[atom1+PADDED_NUM_ATOMS*2], static_cast<unsigned long long>((long long) (force1.z*0x100000000)));
        __threadfence_block();
        atomicAdd(&forceBuffers[atom2],                    static_cast<unsigned long long>((long long) (force2.x*0x100000000)));
        atomicAdd(&forceBuffers[atom2+PADDED_NUM_ATOMS],   static_cast<unsigned long long>((long long) (force2.y*0x100000000)));
        atomicAdd(&forceBuffers[atom2+PADDED_NUM_ATOMS*2], static_cast<unsigned long long>((long long) (force2.z*0x100000000)));
        __threadfence_block();
    }
}


//loop over each angle to calculate chain rule force 
extern "C" __global__ void computeCFluxAngleForce(unsigned long long* __restrict__ forceBuffers, const real4* __restrict__ posq, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, const float3* angleCFluxAngleParams, const float4* angleCFluxBondParams, real* __restrict__ cfPotentials, const uint3* __restrict__ atomIndicesAngle) {


    for (unsigned int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_ANGLES; index += blockDim.x*gridDim.x) {
        uint3 atoms0 = atomIndicesAngle[index];
        unsigned int ia = atoms0.x;
        unsigned int ib = atoms0.y;
        unsigned int ic = atoms0.z;
        real4 pos1 = posq[ia];
        real4 pos2 = posq[ib];
        real4 pos3 = posq[ic];
        
        real3 vba = make_real3(pos1.x-pos2.x, pos1.y-pos2.y, pos1.z-pos2.z); // a-b 
        real3 vbc = make_real3(pos3.x-pos2.x, pos3.y-pos2.y, pos3.z-pos2.z); // c-b
#if APPLY_PERIODIC
        APPLY_PERIODIC_TO_DELTA(vba)
        APPLY_PERIODIC_TO_DELTA(vbc)
#endif
        float3 angleParams = angleCFluxAngleParams[index];
        real jtheta1 = angleParams.y;
        real jtheta2 = angleParams.z;
     
        float4 bondParams = angleCFluxBondParams[index];
        real jbp1  = bondParams.y;
        real jbp2  = bondParams.w;

        real rba2 = vba.x*vba.x + vba.y*vba.y + vba.z*vba.z;
        real rbc2 = vbc.x*vbc.x + vbc.y*vbc.y + vbc.z*vbc.z;

        real rba = SQRT(rba2);
        real rbc = SQRT(rbc2);

        real rba3 = rba2*rba;
        real rbc3 = rbc2*rbc;

        real potA = cfPotentials[ia];
        real potB = cfPotentials[ib];
        real potC = cfPotentials[ic];
         
        real frcxa1 = -(potB-potC)*jbp2*(vba.x)/rba;
        real frcya1 = -(potB-potC)*jbp2*(vba.y)/rba;
        real frcza1 = -(potB-potC)*jbp2*(vba.z)/rba;
        
        real frcxc1 = -(potB-potA)*jbp1*(vbc.x)/rbc;
        real frcyc1 = -(potB-potA)*jbp1*(vbc.y)/rbc;
        real frczc1 = -(potB-potA)*jbp1*(vbc.z)/rbc;
        
        real frcxb1 = -(frcxa1 + frcxc1); 
        real frcyb1 = -(frcya1 + frcyc1);
        real frczb1 = -(frcza1 + frczc1); 

        real dot = vba.x*vbc.x + vba.y*vbc.y + vba.z*vbc.z;
        real term1 = -rba*rbc/SQRT(rba2*rbc2-dot*dot);

        real term2xa = vbc.x/(rba*rbc) - vba.x*dot/(rba3*rbc);
        real term2xc = vba.x/(rba*rbc) - vbc.x*dot/(rba*rbc3);

        real term2ya = vbc.y/(rba*rbc) - vba.y*dot/(rba3*rbc);
        real term2yc = vba.y/(rba*rbc) - vbc.y*dot/(rba*rbc3);

        real term2za = vbc.z/(rba*rbc) - vba.z*dot/(rba3*rbc);
        real term2zc = vba.z/(rba*rbc) - vbc.z*dot/(rba*rbc3);

        real frcxa2 = (potA*jtheta1 - potB*(jtheta1+jtheta2) + potC*jtheta2)*term1*term2xa*RAD_TO_DEG;
        real frcya2 = (potA*jtheta1 - potB*(jtheta1+jtheta2) + potC*jtheta2)*term1*term2ya*RAD_TO_DEG;
        real frcza2 = (potA*jtheta1 - potB*(jtheta1+jtheta2) + potC*jtheta2)*term1*term2za*RAD_TO_DEG;

        real frcxc2 = (potA*jtheta1 - potB*(jtheta1+jtheta2) + potC*jtheta2)*term1*term2xc*RAD_TO_DEG;
        real frcyc2 = (potA*jtheta1 - potB*(jtheta1+jtheta2) + potC*jtheta2)*term1*term2yc*RAD_TO_DEG;
        real frczc2 = (potA*jtheta1 - potB*(jtheta1+jtheta2) + potC*jtheta2)*term1*term2zc*RAD_TO_DEG;

        real frcxb2 = -(frcxa2 + frcxc2);
        real frcyb2 = -(frcya2 + frcyc2);
        real frczb2 = -(frcza2 + frczc2);

        real frcxaT = frcxa1 + frcxa2;
        real frcyaT = frcya1 + frcya2;
        real frczaT = frcza1 + frcza2;
     
        real frcxbT = frcxb1 + frcxb2;
        real frcybT = frcyb1 + frcyb2;
        real frczbT = frczb1 + frczb2;

        real frcxcT = frcxc1 + frcxc2;
        real frcycT = frcyc1 + frcyc2;
        real frczcT = frczc1 + frczc2;

        real3 forceA = make_real3(-frcxaT, -frcyaT, -frczaT);  
        real3 forceB = make_real3(-frcxbT, -frcybT, -frczbT); 
        real3 forceC = make_real3(-frcxcT, -frcycT, -frczcT);

        atomicAdd(&forceBuffers[ia],                    static_cast<unsigned long long>((long long) (forceA.x*0x100000000)));
        atomicAdd(&forceBuffers[ia+PADDED_NUM_ATOMS],   static_cast<unsigned long long>((long long) (forceA.y*0x100000000)));
        atomicAdd(&forceBuffers[ia+PADDED_NUM_ATOMS*2], static_cast<unsigned long long>((long long) (forceA.z*0x100000000)));
        __threadfence_block();
        atomicAdd(&forceBuffers[ib],                    static_cast<unsigned long long>((long long) (forceB.x*0x100000000)));
        atomicAdd(&forceBuffers[ib+PADDED_NUM_ATOMS],   static_cast<unsigned long long>((long long) (forceB.y*0x100000000)));
        atomicAdd(&forceBuffers[ib+PADDED_NUM_ATOMS*2], static_cast<unsigned long long>((long long) (forceB.z*0x100000000)));
        __threadfence_block();
        atomicAdd(&forceBuffers[ic],                    static_cast<unsigned long long>((long long) (forceC.x*0x100000000)));
        atomicAdd(&forceBuffers[ic+PADDED_NUM_ATOMS],   static_cast<unsigned long long>((long long) (forceC.y*0x100000000)));
        atomicAdd(&forceBuffers[ic+PADDED_NUM_ATOMS*2], static_cast<unsigned long long>((long long) (forceC.z*0x100000000)));
        __threadfence_block();
    }
}

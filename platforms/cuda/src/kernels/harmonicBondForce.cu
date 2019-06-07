float2 bondParams = PARAMS[index];
real deltaIdeal = r-bondParams.x;
energy += 0.5f * bondParams.y*deltaIdeal*deltaIdeal;
real dEdR = bondParams.y * deltaIdeal;
real dE = dEdR/r;
real3 ab = make_real3(POSQ[atom1].x - POSQ[atom2].x,POSQ[atom1].y- POSQ[atom2].y, POSQ[atom1].z - POSQ[atom2].z);
#if APPLY_PERIODIC
APPLY_PERIODIC_TO_DELTA(ab);  
#endif
real dEdx = dE * ab.x;
real dEdy = dE * ab.y;
real dEdz = dE * ab.z;
vxx+= ab.x * dEdx;
vxy+= ab.y * dEdx;
vxz+= ab.z * dEdx;
vyy+= ab.y * dEdy;
vyz+= ab.z * dEdy;
vzz+= ab.z * dEdz;

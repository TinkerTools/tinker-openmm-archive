float2 angleParams = PARAMS[index];
real deltaIdeal = theta*RAD_TO_DEG-angleParams.x;
real deltaIdeal2 = deltaIdeal*deltaIdeal;
real deltaIdeal3 = deltaIdeal*deltaIdeal2;
real deltaIdeal4 = deltaIdeal2*deltaIdeal2;
energy += angleParams.y*deltaIdeal2*(1.0f + CUBIC_K*deltaIdeal + QUARTIC_K*deltaIdeal2 + PENTIC_K*deltaIdeal3 + SEXTIC_K*deltaIdeal4);
real dEdAngle = angleParams.y*deltaIdeal*(2.0f + 3.0f*CUBIC_K*deltaIdeal + 4.0f*QUARTIC_K*deltaIdeal2 + 5.0f*PENTIC_K*deltaIdeal3 + 6.0f*SEXTIC_K*deltaIdeal4);
dEdAngle *= RAD_TO_DEG;
#if USES_VIRIAL
float3 ab = make_real3(POSQ[atom1].x - POSQ[atom2].x,POSQ[atom1].y- POSQ[atom2].y, POSQ[atom1].z - POSQ[atom2].z);
float3 cb = make_real3(POSQ[atom3].x - POSQ[atom2].x,POSQ[atom3].y- POSQ[atom2].y, POSQ[atom3].z - POSQ[atom2].z);
#if APPLY_PERIODIC
APPLY_PERIODIC_TO_DELTA(ab);
APPLY_PERIODIC_TO_DELTA(cb);
#endif
float xp= cb.y*ab.z-cb.z*ab.y;
float yp=cb.z*ab.x-cb.x*ab.z;
float zp=cb.x*ab.y-cb.y*ab.x;
float rab2=ab.x*ab.x+ab.y*ab.y+ab.z*ab.z;
float rcb2= cb.x*cb.x+cb.y*cb.y+cb.z*cb.z;
float radp= sqrtf(xp*xp+yp*yp+zp*zp);
float terma = -1.0f*dEdAngle/(rab2*radp);
float termc= dEdAngle/(rcb2*radp);
float dedxia=terma*(ab.y*zp-ab.z*yp);
float dedyia = terma * (ab.z*xp-ab.x*zp);
float dedzia = terma * (ab.x*yp-ab.y*xp);
float dedxic=termc*(cb.y*zp-cb.z*yp);
float dedyic = termc * (cb.z*xp-cb.x*zp);
float dedzic = termc * (cb.x*yp-cb.y*xp);
vxx+=ab.x*dedxia+cb.x*dedxic;
vxy+=ab.y*dedxia+cb.y*dedxic;
vxz+= ab.z*dedxia+cb.z*dedxic;
vyy+=ab.y*dedyia+cb.y*dedyic;
vyz+= ab.z*dedyia + cb.z*dedyic;
vzz += ab.z*dedzia + cb.z*dedzic;
#endif

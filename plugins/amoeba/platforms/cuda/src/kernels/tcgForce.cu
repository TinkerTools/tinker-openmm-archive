#ifdef TCG_POLARIZATION

/**
 * Calculates the polarization energy via -1/2 mu.permField.
 *
 * Parameters:
 * mixed, dimension(*) :: energyBuffer
 * long long, dimension(PADDED_NUM_ATOMS,3) :: permField,permFieldPolar
 * real, dimension(3,*) :: uind,uinp
 */
extern "C" __global__ void tcgPolarizationEnergy(mixed* __restrict__ energyBuffer,
  const long long* __restrict__ permField, const long long* __restrict__ permFieldPolar, const real* __restrict__ uind,
  const real* __restrict__ uinp)
{
  const real fieldScale = 1 / (real)0x100000000;
  int writePos = blockIdx.x * blockDim.x + threadIdx.x;
  for (int index = blockIdx.x * blockDim.x + threadIdx.x; index < NUM_ATOMS; index += blockDim.x * gridDim.x) {
    int idx = 3 * index;
    real edx = permField[index] * fieldScale;
    real edy = permField[index + PADDED_NUM_ATOMS] * fieldScale;
    real edz = permField[index + PADDED_NUM_ATOMS * 2] * fieldScale;
    real epx = permFieldPolar[index] * fieldScale;
    real epy = permFieldPolar[index + PADDED_NUM_ATOMS] * fieldScale;
    real epz = permFieldPolar[index + PADDED_NUM_ATOMS * 2] * fieldScale;
    mixed eterm = -0.25 * ENERGY_SCALE_FACTOR
      * (uind[idx] * epx + uind[idx + 1] * epy + uind[idx + 2] * epz + uinp[idx] * edx + uinp[idx + 1] * edy
          + uinp[idx + 2] * edz);
    energyBuffer[writePos] += eterm;
  }
}

typedef struct {
  real3 pos, force, torque, inducedDipole, inducedDipolePolar, sphericalDipole;
  real q;
  float thole, damp;
#  ifdef INCLUDE_QUADRUPOLES
  real sphericalQuadrupole[5];
#  endif
} TCGPoleData;

inline __device__ void loadTCGPoleData(TCGPoleData& data, int atom, const real4* __restrict__ posq,
  const real* __restrict__ sphericalDipole, const real* __restrict__ sphericalQuadrupole,
  const real* __restrict__ inducedDipole, const real* __restrict__ inducedDipolePolar,
  const float2* __restrict__ dampingAndThole)
{
  real4 atomPosq = posq[atom];
  data.pos = make_real3(atomPosq.x, atomPosq.y, atomPosq.z);
  data.q = atomPosq.w;
  data.sphericalDipole.x = sphericalDipole[atom * 3];
  data.sphericalDipole.y = sphericalDipole[atom * 3 + 1];
  data.sphericalDipole.z = sphericalDipole[atom * 3 + 2];
#  ifdef INCLUDE_QUADRUPOLES
  data.sphericalQuadrupole[0] = sphericalQuadrupole[atom * 5];
  data.sphericalQuadrupole[1] = sphericalQuadrupole[atom * 5 + 1];
  data.sphericalQuadrupole[2] = sphericalQuadrupole[atom * 5 + 2];
  data.sphericalQuadrupole[3] = sphericalQuadrupole[atom * 5 + 3];
  data.sphericalQuadrupole[4] = sphericalQuadrupole[atom * 5 + 4];
#  endif
  data.inducedDipole.x = inducedDipole[atom * 3];
  data.inducedDipole.y = inducedDipole[atom * 3 + 1];
  data.inducedDipole.z = inducedDipole[atom * 3 + 2];
  data.inducedDipolePolar.x = inducedDipolePolar[atom * 3];
  data.inducedDipolePolar.y = inducedDipolePolar[atom * 3 + 1];
  data.inducedDipolePolar.z = inducedDipolePolar[atom * 3 + 2];
  float2 temp = dampingAndThole[atom];
  data.damp = temp.x;
  data.thole = temp.y;
}

inline __device__ real computeDScaleFactor(unsigned int polarizationGroup, int index)
{
  return (polarizationGroup & 1 << index ? 0 : 1);
}

inline __device__ float computeMScaleFactor(uint2 covalent, int index)
{
  int mask = 1 << index;
  bool x = (covalent.x & mask);
  bool y = (covalent.y & mask);
  return (x ? (y ? 0.0f : 0.4f) : (y ? 0.8f : 1.0f));
}

inline __device__ float computePScaleFactor(uint2 covalent, unsigned int polarizationGroup, int index)
{
  int mask = 1 << index;
  bool x = (covalent.x & mask);
  bool y = (covalent.y & mask);
  bool p = (polarizationGroup & mask);
  return (x && y ? 0.0f : (x && p ? 0.5f : 1.0f));
}

// self multipole-multipole energy and self multipole-induced dipole torque
inline __device__ void tcgSelf(TCGPoleData& atom1, mixed& energy)
{
  real cii = atom1.q * atom1.q;
  real dii = dot(atom1.sphericalDipole, atom1.sphericalDipole);
#  ifdef INCLUDE_QUADRUPOLES
  real qii = (atom1.sphericalQuadrupole[0] * atom1.sphericalQuadrupole[0]
    + atom1.sphericalQuadrupole[1] * atom1.sphericalQuadrupole[1]
    + atom1.sphericalQuadrupole[2] * atom1.sphericalQuadrupole[2]
    + atom1.sphericalQuadrupole[3] * atom1.sphericalQuadrupole[3]
    + atom1.sphericalQuadrupole[4] * atom1.sphericalQuadrupole[4]);
#  else
  real qii = 0;
#  endif
  real prefac = -EWALD_ALPHA / SQRT_PI;
  real a2 = EWALD_ALPHA * EWALD_ALPHA;
  real a4 = a2 * a2;
  energy += prefac * (cii + ((real)2 / 3) * a2 * dii + ((real)4 / 15) * a4 * qii);

  // self-torque
  real3 dipole = make_real3(atom1.sphericalDipole.y, atom1.sphericalDipole.z, atom1.sphericalDipole.x);
  real3 ui = atom1.inducedDipole + atom1.inducedDipolePolar;
  atom1.torque += ((2 / (real)3) * (EWALD_ALPHA * EWALD_ALPHA * EWALD_ALPHA) / SQRT_PI) * cross(dipole, ui);
}

// real space and self multipole-multipole force and energy,
// and real space multipole-induced dipole force (no energy)
__device__ void tcgPoleRealSpacePairIxn(TCGPoleData& atom1, TCGPoleData& atom2, bool hasExclusions, float dScale,
  float pScale, float mScale, float forceFactor, mixed& energy, real4 periodicBoxSize, real4 invPeriodicBoxSize,
  real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ)
{
  // Compute the displacement.

  real3 delta;
  delta.x = atom2.pos.x - atom1.pos.x;
  delta.y = atom2.pos.y - atom1.pos.y;
  delta.z = atom2.pos.z - atom1.pos.z;
  APPLY_PERIODIC_TO_DELTA(delta)
  real r2 = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
  if (r2 > CUTOFF_SQUARED)
    return;

  real rInv = RSQRT(r2);
  real r = r2 * rInv;

  // Rotate the various dipoles and quadrupoles.

  real qiRotationMatrix[3][3];
  buildQIRotationMatrix(delta, rInv, qiRotationMatrix);

  real3 rotatedDipole1 = rotateDipole(atom1.sphericalDipole, qiRotationMatrix);
  real3 rotatedDipole2 = rotateDipole(atom2.sphericalDipole, qiRotationMatrix);
  real rotatedQuadrupole1[] = {0, 0, 0, 0, 0};
  real rotatedQuadrupole2[] = {0, 0, 0, 0, 0};
#  ifdef INCLUDE_QUADRUPOLES
  rotateQuadupoles(
    qiRotationMatrix, atom1.sphericalQuadrupole, atom2.sphericalQuadrupole, rotatedQuadrupole1, rotatedQuadrupole2);
#  endif

  // The field derivatives at I due to only permanent moments on J, and vice-versa.
  real Vij[9], Vji[9], VjiR[9], VijR[9];

  // The rInvVec array is defined such that the ith element is R^-i, with the
  // dieleectric constant folded in, to avoid conversions later.
  real rInvVec[7], alphaRVec[8], bVec[5];

  rInvVec[1] = rInv;
  for (int i = 2; i < 7; ++i)
    rInvVec[i] = rInvVec[i - 1] * rInv;

  // The alpharVec array is defined such that the ith element is (alpha R)^i,
  // where kappa (alpha in OpenMM parlance) is the Ewald attenuation parameter.
  real ralpha = EWALD_ALPHA * r;
  real exp2a = EXP(-(ralpha * ralpha));
#  ifdef USE_DOUBLE_PRECISION
  const real erfAlphaR = erf(ralpha);
#  else
  // This approximation for erfc is from Abramowitz and Stegun (1964) p. 299.  They cite the following as
  // the original source: C. Hastings, Jr., Approximations for Digital Computers (1955).  It has a maximum
  // error of 1.5e-7.

  const real t = RECIP(1.0f + 0.3275911f * ralpha);
  const real erfAlphaR = 1
    - (0.254829592f + (-0.284496736f + (1.421413741f + (-1.453152027f + 1.061405429f * t) * t) * t) * t) * t * exp2a;
#  endif
  alphaRVec[1] = ralpha;
  for (int i = 2; i < 8; ++i)
    alphaRVec[i] = alphaRVec[i - 1] * ralpha;
  real X = 2 * exp2a / SQRT_PI;
  int doubleFactorial = 1, facCount = 1;
  real tmp = alphaRVec[1];
  bVec[1] = -erfAlphaR;
  for (int i = 2; i < 5; ++i) {
    bVec[i] = bVec[i - 1] + tmp * X / (real)(doubleFactorial);
    facCount = facCount + 2;
    doubleFactorial = doubleFactorial * facCount;
    tmp *= 2 * alphaRVec[2];
  }

  // Energy

  real ePermCoef, dPermCoef;

  // C-C terms (m=0)
  ePermCoef = rInvVec[1] * (mScale + bVec[2] - alphaRVec[1] * X);
  dPermCoef = -0.5f * (mScale + bVec[2]) * rInvVec[2];
  Vij[0] = ePermCoef * atom2.q;
  Vji[0] = ePermCoef * atom1.q;
  VijR[0] = dPermCoef * atom2.q;
  VjiR[0] = dPermCoef * atom1.q;

  // C-D terms (m=0)
  ePermCoef = rInvVec[2] * (mScale + bVec[2]);
  dPermCoef = -rInvVec[3] * (mScale + bVec[2] + alphaRVec[3] * X);
  Vij[0] += -(ePermCoef * rotatedDipole2.x);
  Vji[1] = -(ePermCoef * atom1.q);
  VijR[0] += -(dPermCoef * rotatedDipole2.x);
  VjiR[1] = -(dPermCoef * atom1.q);
  // D-C terms (m=0)
  Vij[1] = ePermCoef * atom2.q;
  Vji[0] += ePermCoef * rotatedDipole1.x;
  VijR[1] = dPermCoef * atom2.q;
  VjiR[0] += dPermCoef * rotatedDipole1.x;

  // D-D terms (m=0)
  const real twoThirds = (real)2 / 3;
  ePermCoef = -twoThirds * rInvVec[3] * (3 * (mScale + bVec[3]) + alphaRVec[3] * X);
  dPermCoef = rInvVec[4] * (3 * (mScale + bVec[3]) + 2 * alphaRVec[5] * X);
  Vij[1] += ePermCoef * rotatedDipole2.x;
  Vji[1] += ePermCoef * rotatedDipole1.x;
  VijR[1] += dPermCoef * rotatedDipole2.x;
  VjiR[1] += dPermCoef * rotatedDipole1.x;
  // D-D terms (m=1)
  ePermCoef = rInvVec[3] * (mScale + bVec[3] - twoThirds * alphaRVec[3] * X);
  dPermCoef = -1.5f * rInvVec[4] * (mScale + bVec[3]);
  Vij[2] = ePermCoef * rotatedDipole2.y;
  Vji[2] = ePermCoef * rotatedDipole1.y;
  VijR[2] = dPermCoef * rotatedDipole2.y;
  VjiR[2] = dPermCoef * rotatedDipole1.y;
  Vij[3] = ePermCoef * rotatedDipole2.z;
  Vji[3] = ePermCoef * rotatedDipole1.z;
  VijR[3] = dPermCoef * rotatedDipole2.z;
  VjiR[3] = dPermCoef * rotatedDipole1.z;

  // C-Q terms (m=0)
  ePermCoef = (mScale + bVec[3]) * rInvVec[3];
  dPermCoef = -((real)1 / 3) * rInvVec[4] * (4.5f * (mScale + bVec[3]) + 2 * alphaRVec[5] * X);
  Vij[0] += ePermCoef * rotatedQuadrupole2[0];
  Vji[4] = ePermCoef * atom1.q;
  VijR[0] += dPermCoef * rotatedQuadrupole2[0];
  VjiR[4] = dPermCoef * atom1.q;
  // Q-C terms (m=0)
  Vij[4] = ePermCoef * atom2.q;
  Vji[0] += ePermCoef * rotatedQuadrupole1[0];
  VijR[4] = dPermCoef * atom2.q;
  VjiR[0] += dPermCoef * rotatedQuadrupole1[0];

  // D-Q terms (m=0)
  const real fourThirds = (real)4 / 3;
  ePermCoef = rInvVec[4] * (3 * (mScale + bVec[3]) + fourThirds * alphaRVec[5] * X);
  dPermCoef = -fourThirds * rInvVec[5] * (4.5f * (mScale + bVec[3]) + (1 + alphaRVec[2]) * alphaRVec[5] * X);
  Vij[1] += ePermCoef * rotatedQuadrupole2[0];
  Vji[4] += ePermCoef * rotatedDipole1.x;
  VijR[1] += dPermCoef * rotatedQuadrupole2[0];
  VjiR[4] += dPermCoef * rotatedDipole1.x;
  // Q-D terms (m=0)
  Vij[4] += -(ePermCoef * rotatedDipole2.x);
  Vji[1] += -(ePermCoef * rotatedQuadrupole1[0]);
  VijR[4] += -(dPermCoef * rotatedDipole2.x);
  VjiR[1] += -(dPermCoef * rotatedQuadrupole1[0]);

  // D-Q terms (m=1)
  const real sqrtThree = SQRT((real)3);
  ePermCoef = -sqrtThree * rInvVec[4] * (mScale + bVec[3]);
  const real fourSqrtOneThird = 4 / sqrt((real)3);
  dPermCoef = fourSqrtOneThird * rInvVec[5] * (1.5f * (mScale + bVec[3]) + 0.5f * alphaRVec[5] * X);
  Vij[2] += ePermCoef * rotatedQuadrupole2[1];
  Vji[5] = ePermCoef * rotatedDipole1.y;
  VijR[2] += dPermCoef * rotatedQuadrupole2[1];
  VjiR[5] = dPermCoef * rotatedDipole1.y;
  Vij[3] += ePermCoef * rotatedQuadrupole2[2];
  Vji[6] = ePermCoef * rotatedDipole1.z;
  VijR[3] += dPermCoef * rotatedQuadrupole2[2];
  VjiR[6] = dPermCoef * rotatedDipole1.z;
  // D-Q terms (m=1)
  Vij[5] = -(ePermCoef * rotatedDipole2.y);
  Vji[2] += -(ePermCoef * rotatedQuadrupole1[1]);
  VijR[5] = -(dPermCoef * rotatedDipole2.y);
  VjiR[2] += -(dPermCoef * rotatedQuadrupole1[1]);
  Vij[6] = -(ePermCoef * rotatedDipole2.z);
  Vji[3] += -(ePermCoef * rotatedQuadrupole1[2]);
  VijR[6] = -(dPermCoef * rotatedDipole2.z);
  VjiR[3] += -(dPermCoef * rotatedQuadrupole1[2]);

  // Q-Q terms (m=0)
  ePermCoef = rInvVec[5] * (6 * (mScale + bVec[4]) + ((real)4 / 45) * (-3 + 10 * alphaRVec[2]) * alphaRVec[5] * X);
  dPermCoef = -rInvVec[6] * (135 * (mScale + bVec[4]) + 4 * (1 + 2 * alphaRVec[2]) * alphaRVec[7] * X) / 9;
  Vij[4] += ePermCoef * rotatedQuadrupole2[0];
  Vji[4] += ePermCoef * rotatedQuadrupole1[0];
  VijR[4] += dPermCoef * rotatedQuadrupole2[0];
  VjiR[4] += dPermCoef * rotatedQuadrupole1[0];
  // Q-Q terms (m=1)
  const real fourOverFifteen = (real)4 / 15;
  ePermCoef = -fourOverFifteen * rInvVec[5] * (15 * (mScale + bVec[4]) + alphaRVec[5] * X);
  dPermCoef = rInvVec[6] * (10 * (mScale + bVec[4]) + fourThirds * alphaRVec[7] * X);
  Vij[5] += ePermCoef * rotatedQuadrupole2[1];
  Vji[5] += ePermCoef * rotatedQuadrupole1[1];
  VijR[5] += dPermCoef * rotatedQuadrupole2[1];
  VjiR[5] += dPermCoef * rotatedQuadrupole1[1];
  Vij[6] += ePermCoef * rotatedQuadrupole2[2];
  Vji[6] += ePermCoef * rotatedQuadrupole1[2];
  VijR[6] += dPermCoef * rotatedQuadrupole2[2];
  VjiR[6] += dPermCoef * rotatedQuadrupole1[2];
  // Q-Q terms (m=2)
  ePermCoef = rInvVec[5] * (mScale + bVec[4] - fourOverFifteen * alphaRVec[5] * X);
  dPermCoef = -2.5f * (mScale + bVec[4]) * rInvVec[6];
  Vij[7] = ePermCoef * rotatedQuadrupole2[3];
  Vji[7] = ePermCoef * rotatedQuadrupole1[3];
  VijR[7] = dPermCoef * rotatedQuadrupole2[3];
  VjiR[7] = dPermCoef * rotatedQuadrupole1[3];
  Vij[8] = ePermCoef * rotatedQuadrupole2[4];
  Vji[8] = ePermCoef * rotatedQuadrupole1[4];
  VijR[8] = dPermCoef * rotatedQuadrupole2[4];
  VjiR[8] = dPermCoef * rotatedQuadrupole1[4];

  energy += forceFactor * 0.5f
    * (atom1.q * Vij[0] + rotatedDipole1.x * Vij[1] + rotatedDipole1.y * Vij[2] + rotatedDipole1.z * Vij[3]
        + rotatedQuadrupole1[0] * Vij[4] + rotatedQuadrupole1[1] * Vij[5] + rotatedQuadrupole1[2] * Vij[6]
        + rotatedQuadrupole1[3] * Vij[7] + rotatedQuadrupole1[4] * Vij[8] + atom2.q * Vji[0] + rotatedDipole2.x * Vji[1]
        + rotatedDipole2.y * Vji[2] + rotatedDipole2.z * Vji[3] + rotatedQuadrupole2[0] * Vji[4]
        + rotatedQuadrupole2[1] * Vji[5] + rotatedQuadrupole2[2] * Vji[6] + rotatedQuadrupole2[3] * Vji[7]
        + rotatedQuadrupole2[4] * Vji[8]);

  // Force

  real3 qiUindI = 0.5f
    * make_real3(qiRotationMatrix[0][1] * atom1.inducedDipole.x + qiRotationMatrix[0][2] * atom1.inducedDipole.y
          + qiRotationMatrix[0][0] * atom1.inducedDipole.z,
        qiRotationMatrix[1][1] * atom1.inducedDipole.x + qiRotationMatrix[1][2] * atom1.inducedDipole.y
          + qiRotationMatrix[1][0] * atom1.inducedDipole.z,
        qiRotationMatrix[2][1] * atom1.inducedDipole.x + qiRotationMatrix[2][2] * atom1.inducedDipole.y
          + qiRotationMatrix[2][0] * atom1.inducedDipole.z);
  real3 qiUindJ = 0.5f
    * make_real3(qiRotationMatrix[0][1] * atom2.inducedDipole.x + qiRotationMatrix[0][2] * atom2.inducedDipole.y
          + qiRotationMatrix[0][0] * atom2.inducedDipole.z,
        qiRotationMatrix[1][1] * atom2.inducedDipole.x + qiRotationMatrix[1][2] * atom2.inducedDipole.y
          + qiRotationMatrix[1][0] * atom2.inducedDipole.z,
        qiRotationMatrix[2][1] * atom2.inducedDipole.x + qiRotationMatrix[2][2] * atom2.inducedDipole.y
          + qiRotationMatrix[2][0] * atom2.inducedDipole.z);
  real3 qiUinpI = 0.5f
    * make_real3(qiRotationMatrix[0][1] * atom1.inducedDipolePolar.x
          + qiRotationMatrix[0][2] * atom1.inducedDipolePolar.y + qiRotationMatrix[0][0] * atom1.inducedDipolePolar.z,
        qiRotationMatrix[1][1] * atom1.inducedDipolePolar.x + qiRotationMatrix[1][2] * atom1.inducedDipolePolar.y
          + qiRotationMatrix[1][0] * atom1.inducedDipolePolar.z,
        qiRotationMatrix[2][1] * atom1.inducedDipolePolar.x + qiRotationMatrix[2][2] * atom1.inducedDipolePolar.y
          + qiRotationMatrix[2][0] * atom1.inducedDipolePolar.z);
  real3 qiUinpJ = 0.5f
    * make_real3(qiRotationMatrix[0][1] * atom2.inducedDipolePolar.x
          + qiRotationMatrix[0][2] * atom2.inducedDipolePolar.y + qiRotationMatrix[0][0] * atom2.inducedDipolePolar.z,
        qiRotationMatrix[1][1] * atom2.inducedDipolePolar.x + qiRotationMatrix[1][2] * atom2.inducedDipolePolar.y
          + qiRotationMatrix[1][0] * atom2.inducedDipolePolar.z,
        qiRotationMatrix[2][1] * atom2.inducedDipolePolar.x + qiRotationMatrix[2][2] * atom2.inducedDipolePolar.y
          + qiRotationMatrix[2][0] * atom2.inducedDipolePolar.z);

  real Vijp[3], Vijd[3], Vjip[3], Vjid[3];
  real eUIndCoef, dUIndCoef, eUInpCoef, dUInpCoef;

  real dmp = atom1.damp * atom2.damp;
  real a = min(atom1.thole, atom2.thole);
  real u = fabs(dmp) > 1.0e-5f ? r / dmp : 1e10f;
  real au3 = a * u * u * u;
  real expau3 = au3 < 50 ? EXP(-au3) : 0;
  real a2u6 = au3 * au3;
  real a3u9 = a2u6 * au3;
  // Thole damping factors for energies
  real thole_c = 1 - expau3;
  real thole_d0 = 1 - expau3 * (1 + 1.5f * au3);
  real thole_d1 = 1 - expau3;
  real thole_q0 = 1 - expau3 * (1 + au3 + a2u6);
  real thole_q1 = 1 - expau3 * (1 + au3);
  // Thole damping factors for derivatives
  real dthole_c = 1 - expau3 * (1 + 1.5f * au3);
  real dthole_d0 = 1 - expau3 * (1 + au3 + 1.5f * a2u6);
  real dthole_d1 = 1 - expau3 * (1 + au3);
  real dthole_q0 = 1 - expau3 * (1 + au3 + 0.25f * a2u6 + 0.75f * a3u9);
  real dthole_q1 = 1 - expau3 * (1 + au3 + 0.75f * a2u6);

  // C-Uind terms (m=0)
  eUIndCoef = rInvVec[2] * (pScale * thole_c + bVec[2]);
  eUInpCoef = rInvVec[2] * (dScale * thole_c + bVec[2]);
  dUIndCoef = -2 * rInvVec[3] * (pScale * dthole_c + bVec[2] + alphaRVec[3] * X);
  dUInpCoef = -2 * rInvVec[3] * (dScale * dthole_c + bVec[2] + alphaRVec[3] * X);
  Vij[0] += -(eUIndCoef * qiUindJ.x + eUInpCoef * qiUinpJ.x);
  VijR[0] += -(dUIndCoef * qiUindJ.x + dUInpCoef * qiUinpJ.x);
  Vjip[0] = -(eUInpCoef * atom1.q);
  Vjid[0] = -(eUIndCoef * atom1.q);
  // Uind-C terms (m=0)
  Vji[0] += eUIndCoef * qiUindI.x + eUInpCoef * qiUinpI.x;
  VjiR[0] += dUIndCoef * qiUindI.x + dUInpCoef * qiUinpI.x;
  Vijp[0] = eUInpCoef * atom2.q;
  Vijd[0] = eUIndCoef * atom2.q;

  // D-Uind terms (m=0)
  eUIndCoef = -twoThirds * rInvVec[3] * (3 * (pScale * thole_d0 + bVec[3]) + alphaRVec[3] * X);
  eUInpCoef = -twoThirds * rInvVec[3] * (3 * (dScale * thole_d0 + bVec[3]) + alphaRVec[3] * X);
  dUIndCoef = rInvVec[4] * (6 * (pScale * dthole_d0 + bVec[3]) + 4 * alphaRVec[5] * X);
  dUInpCoef = rInvVec[4] * (6 * (dScale * dthole_d0 + bVec[3]) + 4 * alphaRVec[5] * X);
  Vij[1] += eUIndCoef * qiUindJ.x + eUInpCoef * qiUinpJ.x;
  Vji[1] += eUIndCoef * qiUindI.x + eUInpCoef * qiUinpI.x;
  VijR[1] += dUIndCoef * qiUindJ.x + dUInpCoef * qiUinpJ.x;
  VjiR[1] += dUIndCoef * qiUindI.x + dUInpCoef * qiUinpI.x;
  Vijp[0] += eUInpCoef * rotatedDipole2.x;
  Vijd[0] += eUIndCoef * rotatedDipole2.x;
  Vjip[0] += eUInpCoef * rotatedDipole1.x;
  Vjid[0] += eUIndCoef * rotatedDipole1.x;
  // D-Uind terms (m=1)
  eUIndCoef = rInvVec[3] * (pScale * thole_d1 + bVec[3] - twoThirds * alphaRVec[3] * X);
  eUInpCoef = rInvVec[3] * (dScale * thole_d1 + bVec[3] - twoThirds * alphaRVec[3] * X);
  dUIndCoef = -3 * rInvVec[4] * (pScale * dthole_d1 + bVec[3]);
  dUInpCoef = -3 * rInvVec[4] * (dScale * dthole_d1 + bVec[3]);
  Vij[2] += eUIndCoef * qiUindJ.y + eUInpCoef * qiUinpJ.y;
  Vji[2] += eUIndCoef * qiUindI.y + eUInpCoef * qiUinpI.y;
  VijR[2] += dUIndCoef * qiUindJ.y + dUInpCoef * qiUinpJ.y;
  VjiR[2] += dUIndCoef * qiUindI.y + dUInpCoef * qiUinpI.y;
  Vij[3] += eUIndCoef * qiUindJ.z + eUInpCoef * qiUinpJ.z;
  Vji[3] += eUIndCoef * qiUindI.z + eUInpCoef * qiUinpI.z;
  VijR[3] += dUIndCoef * qiUindJ.z + dUInpCoef * qiUinpJ.z;
  VjiR[3] += dUIndCoef * qiUindI.z + dUInpCoef * qiUinpI.z;
  Vijp[1] = eUInpCoef * rotatedDipole2.y;
  Vijd[1] = eUIndCoef * rotatedDipole2.y;
  Vjip[1] = eUInpCoef * rotatedDipole1.y;
  Vjid[1] = eUIndCoef * rotatedDipole1.y;
  Vijp[2] = eUInpCoef * rotatedDipole2.z;
  Vijd[2] = eUIndCoef * rotatedDipole2.z;
  Vjip[2] = eUInpCoef * rotatedDipole1.z;
  Vjid[2] = eUIndCoef * rotatedDipole1.z;

  // Uind-Q terms (m=0)
  eUIndCoef = rInvVec[4] * (3 * (pScale * thole_q0 + bVec[3]) + fourThirds * alphaRVec[5] * X);
  eUInpCoef = rInvVec[4] * (3 * (dScale * thole_q0 + bVec[3]) + fourThirds * alphaRVec[5] * X);
  dUIndCoef
    = -fourThirds * rInvVec[5] * (9 * (pScale * dthole_q0 + bVec[3]) + 2 * (1 + alphaRVec[2]) * alphaRVec[5] * X);
  dUInpCoef
    = -fourThirds * rInvVec[5] * (9 * (dScale * dthole_q0 + bVec[3]) + 2 * (1 + alphaRVec[2]) * alphaRVec[5] * X);
  Vji[4] += eUIndCoef * qiUindI.x + eUInpCoef * qiUinpI.x;
  VjiR[4] += dUIndCoef * qiUindI.x + dUInpCoef * qiUinpI.x;
  Vijp[0] += eUInpCoef * rotatedQuadrupole2[0];
  Vijd[0] += eUIndCoef * rotatedQuadrupole2[0];
  // Q-Uind terms (m=0)
  Vij[4] += -(eUIndCoef * qiUindJ.x + eUInpCoef * qiUinpJ.x);
  VijR[4] += -(dUIndCoef * qiUindJ.x + dUInpCoef * qiUinpJ.x);
  Vjip[0] += -(eUInpCoef * rotatedQuadrupole1[0]);
  Vjid[0] += -(eUIndCoef * rotatedQuadrupole1[0]);

  // Uind-Q terms (m=1)
  eUIndCoef = -sqrtThree * rInvVec[4] * (pScale * thole_q1 + bVec[3]);
  eUInpCoef = -sqrtThree * rInvVec[4] * (dScale * thole_q1 + bVec[3]);
  dUIndCoef = fourSqrtOneThird * rInvVec[5] * (3 * (pScale * dthole_q1 + bVec[3]) + alphaRVec[5] * X);
  dUInpCoef = fourSqrtOneThird * rInvVec[5] * (3 * (dScale * dthole_q1 + bVec[3]) + alphaRVec[5] * X);
  Vji[5] += eUIndCoef * qiUindI.y + eUInpCoef * qiUinpI.y;
  VjiR[5] += dUIndCoef * qiUindI.y + dUInpCoef * qiUinpI.y;
  Vji[6] += eUIndCoef * qiUindI.z + eUInpCoef * qiUinpI.z;
  VjiR[6] += dUIndCoef * qiUindI.z + dUInpCoef * qiUinpI.z;
  Vijp[1] += eUInpCoef * rotatedQuadrupole2[1];
  Vijd[1] += eUIndCoef * rotatedQuadrupole2[1];
  Vijp[2] += eUInpCoef * rotatedQuadrupole2[2];
  Vijd[2] += eUIndCoef * rotatedQuadrupole2[2];
  // Q-Uind terms (m=1)
  Vij[5] += -(eUIndCoef * qiUindJ.y + eUInpCoef * qiUinpJ.y);
  VijR[5] += -(dUIndCoef * qiUindJ.y + dUInpCoef * qiUinpJ.y);
  Vij[6] += -(eUIndCoef * qiUindJ.z + eUInpCoef * qiUinpJ.z);
  VijR[6] += -(dUIndCoef * qiUindJ.z + dUInpCoef * qiUinpJ.z);
  Vjip[1] += -(eUInpCoef * rotatedQuadrupole1[1]);
  Vjid[1] += -(eUIndCoef * rotatedQuadrupole1[1]);
  Vjip[2] += -(eUInpCoef * rotatedQuadrupole1[2]);
  Vjid[2] += -(eUIndCoef * rotatedQuadrupole1[2]);

  real fIZ = atom1.q * VijR[0] + rotatedDipole1.x * VijR[1] + rotatedDipole1.y * VijR[2] + rotatedDipole1.z * VijR[3]
    + rotatedQuadrupole1[0] * VijR[4] + rotatedQuadrupole1[1] * VijR[5] + rotatedQuadrupole1[2] * VijR[6]
    + rotatedQuadrupole1[3] * VijR[7] + rotatedQuadrupole1[4] * VijR[8];
  real fJZ = atom2.q * VjiR[0] + rotatedDipole2.x * VjiR[1] + rotatedDipole2.y * VjiR[2] + rotatedDipole2.z * VjiR[3]
    + rotatedQuadrupole2[0] * VjiR[4] + rotatedQuadrupole2[1] * VjiR[5] + rotatedQuadrupole2[2] * VjiR[6]
    + rotatedQuadrupole2[3] * VjiR[7] + rotatedQuadrupole2[4] * VjiR[8];
  real EIX = rotatedDipole1.z * Vij[1] - rotatedDipole1.x * Vij[3] + sqrtThree * rotatedQuadrupole1[2] * Vij[4]
    + rotatedQuadrupole1[4] * Vij[5] - (sqrtThree * rotatedQuadrupole1[0] + rotatedQuadrupole1[3]) * Vij[6]
    + rotatedQuadrupole1[2] * Vij[7] - rotatedQuadrupole1[1] * Vij[8];
  real EIY = -rotatedDipole1.y * Vij[1] + rotatedDipole1.x * Vij[2] - sqrtThree * rotatedQuadrupole1[1] * Vij[4]
    + (sqrtThree * rotatedQuadrupole1[0] - rotatedQuadrupole1[3]) * Vij[5] - rotatedQuadrupole1[4] * Vij[6]
    + rotatedQuadrupole1[1] * Vij[7] + rotatedQuadrupole1[2] * Vij[8];
  real EIZ = -rotatedDipole1.z * Vij[2] + rotatedDipole1.y * Vij[3] - rotatedQuadrupole1[2] * Vij[5]
    + rotatedQuadrupole1[1] * Vij[6] - 2 * rotatedQuadrupole1[4] * Vij[7] + 2 * rotatedQuadrupole1[3] * Vij[8];
  real EJX = rotatedDipole2.z * Vji[1] - rotatedDipole2.x * Vji[3] + sqrtThree * rotatedQuadrupole2[2] * Vji[4]
    + rotatedQuadrupole2[4] * Vji[5] - (sqrtThree * rotatedQuadrupole2[0] + rotatedQuadrupole2[3]) * Vji[6]
    + rotatedQuadrupole2[2] * Vji[7] - rotatedQuadrupole2[1] * Vji[8];
  real EJY = -rotatedDipole2.y * Vji[1] + rotatedDipole2.x * Vji[2] - sqrtThree * rotatedQuadrupole2[1] * Vji[4]
    + (sqrtThree * rotatedQuadrupole2[0] - rotatedQuadrupole2[3]) * Vji[5] - rotatedQuadrupole2[4] * Vji[6]
    + rotatedQuadrupole2[1] * Vji[7] + rotatedQuadrupole2[2] * Vji[8];
  real EJZ = -rotatedDipole2.z * Vji[2] + rotatedDipole2.y * Vji[3] - rotatedQuadrupole2[2] * Vji[5]
    + rotatedQuadrupole2[1] * Vji[6] - 2 * rotatedQuadrupole2[4] * Vji[7] + 2 * rotatedQuadrupole2[3] * Vji[8];

  // Define the torque intermediates for the induced dipoles. These are simply the induced dipole torque
  // intermediates dotted with the field due to permanent moments only, at each center. We inline the
  // induced dipole torque intermediates here, for simplicity. N.B. There are no torques on the dipoles
  // themselves, so we accumulate the torque intermediates into separate variables to allow them to be
  // used only in the force calculation.
  //
  // The torque about the x axis (needed to obtain the y force on the induced dipoles, below)
  //    qiUindIx[0] = qiQUindI[2];    qiUindIx[1] = 0;    qiUindIx[2] = -qiQUindI[0]
  real iEIX = qiUinpI.z * Vijp[0] + qiUindI.z * Vijd[0] - qiUinpI.x * Vijp[2] - qiUindI.x * Vijd[2];
  real iEJX = qiUinpJ.z * Vjip[0] + qiUindJ.z * Vjid[0] - qiUinpJ.x * Vjip[2] - qiUindJ.x * Vjid[2];
  // The torque about the y axis (needed to obtain the x force on the induced dipoles, below)
  //    qiUindIy[0] = -qiQUindI[1];   qiUindIy[1] = qiQUindI[0];    qiUindIy[2] = 0
  real iEIY = qiUinpI.x * Vijp[1] + qiUindI.x * Vijd[1] - qiUinpI.y * Vijp[0] - qiUindI.y * Vijd[0];
  real iEJY = qiUinpJ.x * Vjip[1] + qiUindJ.x * Vjid[1] - qiUinpJ.y * Vjip[0] - qiUindJ.y * Vjid[0];

  // The quasi-internal frame forces and torques.  Note that the induced torque intermediates are
  // used in the force expression, but not in the torques; the induced dipoles are isotropic.
  real qiForce[3] = {rInv * (EIY + EJY + iEIY + iEJY), -rInv * (EIX + EJX + iEIX + iEJX), -(fJZ + fIZ)};
  real qiTorqueI[3] = {-EIX, -EIY, -EIZ};
  real qiTorqueJ[3] = {-EJX, -EJY, -EJZ};

  real3 force = make_real3(
    qiRotationMatrix[1][1] * qiForce[0] + qiRotationMatrix[2][1] * qiForce[1] + qiRotationMatrix[0][1] * qiForce[2],
    qiRotationMatrix[1][2] * qiForce[0] + qiRotationMatrix[2][2] * qiForce[1] + qiRotationMatrix[0][2] * qiForce[2],
    qiRotationMatrix[1][0] * qiForce[0] + qiRotationMatrix[2][0] * qiForce[1] + qiRotationMatrix[0][0] * qiForce[2]);
  atom1.force += force;
  atom1.torque += make_real3(qiRotationMatrix[1][1] * qiTorqueI[0] + qiRotationMatrix[2][1] * qiTorqueI[1]
      + qiRotationMatrix[0][1] * qiTorqueI[2],
    qiRotationMatrix[1][2] * qiTorqueI[0] + qiRotationMatrix[2][2] * qiTorqueI[1]
      + qiRotationMatrix[0][2] * qiTorqueI[2],
    qiRotationMatrix[1][0] * qiTorqueI[0] + qiRotationMatrix[2][0] * qiTorqueI[1]
      + qiRotationMatrix[0][0] * qiTorqueI[2]);
  if (forceFactor == 1) {
    atom2.force -= force;
    atom2.torque += make_real3(qiRotationMatrix[1][1] * qiTorqueJ[0] + qiRotationMatrix[2][1] * qiTorqueJ[1]
        + qiRotationMatrix[0][1] * qiTorqueJ[2],
      qiRotationMatrix[1][2] * qiTorqueJ[0] + qiRotationMatrix[2][2] * qiTorqueJ[1]
        + qiRotationMatrix[0][2] * qiTorqueJ[2],
      qiRotationMatrix[1][0] * qiTorqueJ[0] + qiRotationMatrix[2][0] * qiTorqueJ[1]
        + qiRotationMatrix[0][0] * qiTorqueJ[2]);
  }
}

// See also: void computeElectrostatics(...) in pmeMultipoleElectrostatics.cu
// real space and self multipole-multipole force and energy,
// and real space multipole-induced dipole force (no energy)
extern "C" __global__ void tcgSelfAndRealSpace(unsigned long long* __restrict__ forceBuffers,
  unsigned long long* __restrict__ torqueBuffers, mixed* __restrict__ energyBuffer, const real4* __restrict__ posq,
  const uint2* __restrict__ covalentFlags, const unsigned int* __restrict__ polarizationGroupFlags,
  const ushort2* __restrict__ exclusionTiles, unsigned int startTileIndex, unsigned int numTileIndices,
#  ifdef USE_CUTOFF
  const int* __restrict__ tiles, const unsigned int* __restrict__ interactionCount, real4 periodicBoxSize,
  real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles,
  const real4* __restrict__ blockCenter, const unsigned int* __restrict__ interactingAtoms,
#  endif
  const real* __restrict__ sphericalDipole, const real* __restrict__ sphericalQuadrupole,
  const real* __restrict__ inducedDipole, const real* __restrict__ inducedDipolePolar,
  const float2* __restrict__ dampingAndThole)
{
  const unsigned int totalWarps = (blockDim.x * gridDim.x) / TILE_SIZE;
  const unsigned int warp = (blockIdx.x * blockDim.x + threadIdx.x) / TILE_SIZE;
  const unsigned int tgx = threadIdx.x & (TILE_SIZE - 1);
  const unsigned int tbx = threadIdx.x - tgx;
  mixed energy = 0;
  __shared__ TCGPoleData localData[MPOLE_THREAD_BLOCK_SIZE];
  __shared__ int atomIndices[MPOLE_THREAD_BLOCK_SIZE];
  __shared__ volatile int skipTiles[MPOLE_THREAD_BLOCK_SIZE];

  // First loop: process tiles that contain exclusions.

  const unsigned int firstExclusionTile
    = FIRST_EXCLUSION_TILE + warp * (LAST_EXCLUSION_TILE - FIRST_EXCLUSION_TILE) / totalWarps;
  const unsigned int lastExclusionTile
    = FIRST_EXCLUSION_TILE + (warp + 1) * (LAST_EXCLUSION_TILE - FIRST_EXCLUSION_TILE) / totalWarps;
  for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
    const ushort2 tileIndices = exclusionTiles[pos];
    const unsigned int x = tileIndices.x;
    const unsigned int y = tileIndices.y;
    TCGPoleData data;
    unsigned int atom1 = x * TILE_SIZE + tgx;
    loadTCGPoleData(
      data, atom1, posq, sphericalDipole, sphericalQuadrupole, inducedDipole, inducedDipolePolar, dampingAndThole);
    data.force = make_real3(0);
    data.torque = make_real3(0);
    uint2 covalent = covalentFlags[pos * TILE_SIZE + tgx];
    unsigned int polarizationGroup = polarizationGroupFlags[pos * TILE_SIZE + tgx];

    if (x == y) {
      // This tile is on the diagonal.

      localData[threadIdx.x].pos = data.pos;
      localData[threadIdx.x].q = data.q;
      localData[threadIdx.x].sphericalDipole = data.sphericalDipole;
#  ifdef INCLUDE_QUADRUPOLES
      localData[threadIdx.x].sphericalQuadrupole[0] = data.sphericalQuadrupole[0];
      localData[threadIdx.x].sphericalQuadrupole[1] = data.sphericalQuadrupole[1];
      localData[threadIdx.x].sphericalQuadrupole[2] = data.sphericalQuadrupole[2];
      localData[threadIdx.x].sphericalQuadrupole[3] = data.sphericalQuadrupole[3];
      localData[threadIdx.x].sphericalQuadrupole[4] = data.sphericalQuadrupole[4];
#  endif
      localData[threadIdx.x].inducedDipole = data.inducedDipole;
      localData[threadIdx.x].inducedDipolePolar = data.inducedDipolePolar;
      localData[threadIdx.x].thole = data.thole;
      localData[threadIdx.x].damp = data.damp;

      // Compute forces.

      for (unsigned int j = 0; j < TILE_SIZE; j++) {
        int atom2 = y * TILE_SIZE + j;
        if (atom1 != atom2 && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
          float d = computeDScaleFactor(polarizationGroup, j);
          float p = computePScaleFactor(covalent, polarizationGroup, j);
          float m = computeMScaleFactor(covalent, j);
          tcgPoleRealSpacePairIxn(data, localData[tbx + j], true, d, p, m, 0.5f, energy, periodicBoxSize,
            invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);
        }
      }
      if (atom1 < NUM_ATOMS)
        tcgSelf(data, energy);
      data.force *= -ENERGY_SCALE_FACTOR;
      data.torque *= ENERGY_SCALE_FACTOR;
      atomicAdd(&forceBuffers[atom1], static_cast<unsigned long long>((long long)(data.force.x * 0x100000000)));
      atomicAdd(&forceBuffers[atom1 + PADDED_NUM_ATOMS],
        static_cast<unsigned long long>((long long)(data.force.y * 0x100000000)));
      atomicAdd(&forceBuffers[atom1 + 2 * PADDED_NUM_ATOMS],
        static_cast<unsigned long long>((long long)(data.force.z * 0x100000000)));
      atomicAdd(&torqueBuffers[atom1], static_cast<unsigned long long>((long long)(data.torque.x * 0x100000000)));
      atomicAdd(&torqueBuffers[atom1 + PADDED_NUM_ATOMS],
        static_cast<unsigned long long>((long long)(data.torque.y * 0x100000000)));
      atomicAdd(&torqueBuffers[atom1 + 2 * PADDED_NUM_ATOMS],
        static_cast<unsigned long long>((long long)(data.torque.z * 0x100000000)));
    } else {
      // This is an off-diagonal tile.

      unsigned int j = y * TILE_SIZE + tgx;
      loadTCGPoleData(localData[threadIdx.x], j, posq, sphericalDipole, sphericalQuadrupole, inducedDipole,
        inducedDipolePolar, dampingAndThole);
      localData[threadIdx.x].force = make_real3(0);
      localData[threadIdx.x].torque = make_real3(0);
      unsigned int tj = tgx;
      for (j = 0; j < TILE_SIZE; j++) {
        int atom2 = y * TILE_SIZE + tj;
        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
          float d = computeDScaleFactor(polarizationGroup, tj);
          float p = computePScaleFactor(covalent, polarizationGroup, tj);
          float m = computeMScaleFactor(covalent, tj);
          tcgPoleRealSpacePairIxn(data, localData[tbx + tj], true, d, p, m, 1, energy, periodicBoxSize,
            invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);
        }
        tj = (tj + 1) & (TILE_SIZE - 1);
      }
      data.force *= -ENERGY_SCALE_FACTOR;
      data.torque *= ENERGY_SCALE_FACTOR;
      localData[threadIdx.x].force *= -ENERGY_SCALE_FACTOR;
      localData[threadIdx.x].torque *= ENERGY_SCALE_FACTOR;
      unsigned int offset = x * TILE_SIZE + tgx;
      atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long)(data.force.x * 0x100000000)));
      atomicAdd(&forceBuffers[offset + PADDED_NUM_ATOMS],
        static_cast<unsigned long long>((long long)(data.force.y * 0x100000000)));
      atomicAdd(&forceBuffers[offset + 2 * PADDED_NUM_ATOMS],
        static_cast<unsigned long long>((long long)(data.force.z * 0x100000000)));
      atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long)(data.torque.x * 0x100000000)));
      atomicAdd(&torqueBuffers[offset + PADDED_NUM_ATOMS],
        static_cast<unsigned long long>((long long)(data.torque.y * 0x100000000)));
      atomicAdd(&torqueBuffers[offset + 2 * PADDED_NUM_ATOMS],
        static_cast<unsigned long long>((long long)(data.torque.z * 0x100000000)));
      offset = y * TILE_SIZE + tgx;
      atomicAdd(&forceBuffers[offset],
        static_cast<unsigned long long>((long long)(localData[threadIdx.x].force.x * 0x100000000)));
      atomicAdd(&forceBuffers[offset + PADDED_NUM_ATOMS],
        static_cast<unsigned long long>((long long)(localData[threadIdx.x].force.y * 0x100000000)));
      atomicAdd(&forceBuffers[offset + 2 * PADDED_NUM_ATOMS],
        static_cast<unsigned long long>((long long)(localData[threadIdx.x].force.z * 0x100000000)));
      atomicAdd(&torqueBuffers[offset],
        static_cast<unsigned long long>((long long)(localData[threadIdx.x].torque.x * 0x100000000)));
      atomicAdd(&torqueBuffers[offset + PADDED_NUM_ATOMS],
        static_cast<unsigned long long>((long long)(localData[threadIdx.x].torque.y * 0x100000000)));
      atomicAdd(&torqueBuffers[offset + 2 * PADDED_NUM_ATOMS],
        static_cast<unsigned long long>((long long)(localData[threadIdx.x].torque.z * 0x100000000)));
    }
  }

  // Second loop: tiles without exclusions, either from the neighbor list (with cutoff) or just enumerating all
  // of them (no cutoff).

#  ifdef USE_CUTOFF
  const unsigned int numTiles = interactionCount[0];
  if (numTiles > maxTiles)
    return; // There wasn't enough memory for the neighbor list.
  int pos = (int)(numTiles > maxTiles ? startTileIndex + warp * (long long)numTileIndices / totalWarps
                                      : warp * (long long)numTiles / totalWarps);
  int end = (int)(numTiles > maxTiles ? startTileIndex + (warp + 1) * (long long)numTileIndices / totalWarps
                                      : (warp + 1) * (long long)numTiles / totalWarps);
#  else
  const unsigned int numTiles = numTileIndices;
  int pos = (int)(startTileIndex + warp * (long long)numTiles / totalWarps);
  int end = (int)(startTileIndex + (warp + 1) * (long long)numTiles / totalWarps);
#  endif
  int skipBase = 0;
  int currentSkipIndex = tbx;
  skipTiles[threadIdx.x] = -1;

  while (pos < end) {
    bool includeTile = true;

    // Extract the coordinates of this tile.

    int x, y;
#  ifdef USE_CUTOFF
    x = tiles[pos];
#  else
    y = (int)floor(NUM_BLOCKS + 0.5f - SQRT((NUM_BLOCKS + 0.5f) * (NUM_BLOCKS + 0.5f) - 2 * pos));
    x = (pos - y * NUM_BLOCKS + y * (y + 1) / 2);
    if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
      y += (x < y ? -1 : 1);
      x = (pos - y * NUM_BLOCKS + y * (y + 1) / 2);
    }

    // Skip over tiles that have exclusions, since they were already processed.

    while (skipTiles[tbx + TILE_SIZE - 1] < pos) {
      if (skipBase + tgx < NUM_TILES_WITH_EXCLUSIONS) {
        ushort2 tile = exclusionTiles[skipBase + tgx];
        skipTiles[threadIdx.x] = tile.x + tile.y * NUM_BLOCKS - tile.y * (tile.y + 1) / 2;
      } else
        skipTiles[threadIdx.x] = end;
      skipBase += TILE_SIZE;
      currentSkipIndex = tbx;
    }
    while (skipTiles[currentSkipIndex] < pos)
      currentSkipIndex++;
    includeTile = (skipTiles[currentSkipIndex] != pos);
#  endif
    if (includeTile) {
      unsigned int atom1 = x * TILE_SIZE + tgx;

      // Load atom data for this tile.

      TCGPoleData data;
      loadTCGPoleData(
        data, atom1, posq, sphericalDipole, sphericalQuadrupole, inducedDipole, inducedDipolePolar, dampingAndThole);
      data.force = make_real3(0);
      data.torque = make_real3(0);
#  ifdef USE_CUTOFF
      unsigned int j = interactingAtoms[pos * TILE_SIZE + tgx];
#  else
      unsigned int j = y * TILE_SIZE + tgx;
#  endif
      atomIndices[threadIdx.x] = j;
      loadTCGPoleData(localData[threadIdx.x], j, posq, sphericalDipole, sphericalQuadrupole, inducedDipole,
        inducedDipolePolar, dampingAndThole);
      localData[threadIdx.x].force = make_real3(0);
      localData[threadIdx.x].torque = make_real3(0);

      // Compute forces.

      unsigned int tj = tgx;
      for (j = 0; j < TILE_SIZE; j++) {
        int atom2 = atomIndices[tbx + tj];
        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
          tcgPoleRealSpacePairIxn(data, localData[tbx + tj], false, 1, 1, 1, 1, energy, periodicBoxSize,
            invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);
        }
        tj = (tj + 1) & (TILE_SIZE - 1);
      }
      data.force *= -ENERGY_SCALE_FACTOR;
      data.torque *= ENERGY_SCALE_FACTOR;
      localData[threadIdx.x].force *= -ENERGY_SCALE_FACTOR;
      localData[threadIdx.x].torque *= ENERGY_SCALE_FACTOR;

      // Write results.

      unsigned int offset = x * TILE_SIZE + tgx;
      atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long)(data.force.x * 0x100000000)));
      atomicAdd(&forceBuffers[offset + PADDED_NUM_ATOMS],
        static_cast<unsigned long long>((long long)(data.force.y * 0x100000000)));
      atomicAdd(&forceBuffers[offset + 2 * PADDED_NUM_ATOMS],
        static_cast<unsigned long long>((long long)(data.force.z * 0x100000000)));
      atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long)(data.torque.x * 0x100000000)));
      atomicAdd(&torqueBuffers[offset + PADDED_NUM_ATOMS],
        static_cast<unsigned long long>((long long)(data.torque.y * 0x100000000)));
      atomicAdd(&torqueBuffers[offset + 2 * PADDED_NUM_ATOMS],
        static_cast<unsigned long long>((long long)(data.torque.z * 0x100000000)));
#  ifdef USE_CUTOFF
      offset = atomIndices[threadIdx.x];
#  else
      offset = y * TILE_SIZE + tgx;
#  endif
      atomicAdd(&forceBuffers[offset],
        static_cast<unsigned long long>((long long)(localData[threadIdx.x].force.x * 0x100000000)));
      atomicAdd(&forceBuffers[offset + PADDED_NUM_ATOMS],
        static_cast<unsigned long long>((long long)(localData[threadIdx.x].force.y * 0x100000000)));
      atomicAdd(&forceBuffers[offset + 2 * PADDED_NUM_ATOMS],
        static_cast<unsigned long long>((long long)(localData[threadIdx.x].force.z * 0x100000000)));
      atomicAdd(&torqueBuffers[offset],
        static_cast<unsigned long long>((long long)(localData[threadIdx.x].torque.x * 0x100000000)));
      atomicAdd(&torqueBuffers[offset + PADDED_NUM_ATOMS],
        static_cast<unsigned long long>((long long)(localData[threadIdx.x].torque.y * 0x100000000)));
      atomicAdd(&torqueBuffers[offset + 2 * PADDED_NUM_ATOMS],
        static_cast<unsigned long long>((long long)(localData[threadIdx.x].torque.z * 0x100000000)));
    }
    pos++;
  }
  energyBuffer[blockIdx.x * blockDim.x + threadIdx.x] += energy * ENERGY_SCALE_FACTOR;
}

////////////////////////////////////////////////////////////////////////////////

typedef struct {
  real3 pos;
  real3 field, fieldPolar, inducedDipole, inducedDipolePolar;
  // IF CALC GRADIENT
  // real fieldGradient[6], fieldGradientPolar[6];
  // END IF
  float thole, damp;
} AtomData0;

typedef struct {
  real3 pos;
  real3 field, fieldPolar, inducedDipole, inducedDipolePolar;
  // IF CALC GRADIENT
  real fieldGradient[6], fieldGradientPolar[6];
  // END IF
  float thole, damp;
} AtomData1;

inline __device__ void loadAtomData0(AtomData0& data, int atom, const real4* __restrict__ posq,
  const real* __restrict__ inducedDipole, const real* __restrict__ inducedDipolePolar,
  const float2* __restrict__ dampingAndThole)
{
  real4 atomPosq = posq[atom];
  data.pos = make_real3(atomPosq.x, atomPosq.y, atomPosq.z);
  data.inducedDipole.x = inducedDipole[atom * 3];
  data.inducedDipole.y = inducedDipole[atom * 3 + 1];
  data.inducedDipole.z = inducedDipole[atom * 3 + 2];
  data.inducedDipolePolar.x = inducedDipolePolar[atom * 3];
  data.inducedDipolePolar.y = inducedDipolePolar[atom * 3 + 1];
  data.inducedDipolePolar.z = inducedDipolePolar[atom * 3 + 2];
  float2 temp = dampingAndThole[atom];
  data.damp = temp.x;
  data.thole = temp.y;
}

inline __device__ void loadAtomData1(AtomData1& data, int atom, const real4* __restrict__ posq,
  const real* __restrict__ inducedDipole, const real* __restrict__ inducedDipolePolar,
  const float2* __restrict__ dampingAndThole)
{
  real4 atomPosq = posq[atom];
  data.pos = make_real3(atomPosq.x, atomPosq.y, atomPosq.z);
  data.inducedDipole.x = inducedDipole[atom * 3];
  data.inducedDipole.y = inducedDipole[atom * 3 + 1];
  data.inducedDipole.z = inducedDipole[atom * 3 + 2];
  data.inducedDipolePolar.x = inducedDipolePolar[atom * 3];
  data.inducedDipolePolar.y = inducedDipolePolar[atom * 3 + 1];
  data.inducedDipolePolar.z = inducedDipolePolar[atom * 3 + 2];
  float2 temp = dampingAndThole[atom];
  data.damp = temp.x;
  data.thole = temp.y;
}

inline __device__ void zeroAtomData0(AtomData0& data)
{
  data.field = make_real3(0);
  data.fieldPolar = make_real3(0);
}

inline __device__ void zeroAtomData1(AtomData1& data)
{
  data.field = make_real3(0);
  data.fieldPolar = make_real3(0);

  for (int i = 0; i < 6; ++i) {
    data.fieldGradient[i] = 0;
    data.fieldGradientPolar[i] = 0;
  }
}

inline __device__ void saveAtomData0(
  int index, AtomData0& data, unsigned long long* __restrict__ field, unsigned long long* __restrict__ fieldPolar)
{
  atomicAdd(&field[index], static_cast<unsigned long long>((long long)(data.field.x * 0x100000000)));
  atomicAdd(&field[index + PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long)(data.field.y * 0x100000000)));
  atomicAdd(
    &field[index + 2 * PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long)(data.field.z * 0x100000000)));
  atomicAdd(&fieldPolar[index], static_cast<unsigned long long>((long long)(data.fieldPolar.x * 0x100000000)));
  atomicAdd(&fieldPolar[index + PADDED_NUM_ATOMS],
    static_cast<unsigned long long>((long long)(data.fieldPolar.y * 0x100000000)));
  atomicAdd(&fieldPolar[index + 2 * PADDED_NUM_ATOMS],
    static_cast<unsigned long long>((long long)(data.fieldPolar.z * 0x100000000)));
}

inline __device__ void saveAtomData1(int index, AtomData1& data, unsigned long long* __restrict__ field,
  unsigned long long* __restrict__ fieldPolar, unsigned long long* __restrict__ fieldGradient,
  unsigned long long* __restrict__ fieldGradientPolar)
{
  atomicAdd(&field[index], static_cast<unsigned long long>((long long)(data.field.x * 0x100000000)));
  atomicAdd(&field[index + PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long)(data.field.y * 0x100000000)));
  atomicAdd(
    &field[index + 2 * PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long)(data.field.z * 0x100000000)));
  atomicAdd(&fieldPolar[index], static_cast<unsigned long long>((long long)(data.fieldPolar.x * 0x100000000)));
  atomicAdd(&fieldPolar[index + PADDED_NUM_ATOMS],
    static_cast<unsigned long long>((long long)(data.fieldPolar.y * 0x100000000)));
  atomicAdd(&fieldPolar[index + 2 * PADDED_NUM_ATOMS],
    static_cast<unsigned long long>((long long)(data.fieldPolar.z * 0x100000000)));

  for (int i = 0; i < 6; i++) {
    atomicAdd(
      &fieldGradient[6 * index + i], static_cast<unsigned long long>((long long)(data.fieldGradient[i] * 0x100000000)));
    atomicAdd(&fieldGradientPolar[6 * index + i],
      static_cast<unsigned long long>((long long)(data.fieldGradientPolar[i] * 0x100000000)));
  }
}

#  define SAVE_ATOM_DATA0(index, data) saveAtomData0(index, data, field, fieldPolar);

#  define SAVE_ATOM_DATA1(index, data) saveAtomData1(index, data, field, fieldPolar, fieldGradient, fieldGradientPolar);

#  ifdef USE_EWALD
// see also: void computeOneInteraction() in multipoleInducedField.cu
__device__ void tcgPairField0(AtomData0& atom1, AtomData0& atom2, real3 deltaR, bool isSelfInteraction)
{
  if (isSelfInteraction)
    return;
  real scale1, scale2;
  real r2 = dot(deltaR, deltaR);
  if (r2 < CUTOFF_SQUARED) {
    real rI = RSQRT(r2);
    real r = RECIP(rI);
    real rI2 = rI * rI;

    // calculate the error function damping terms

    real ralpha = EWALD_ALPHA * r;
    real exp2a = EXP(-(ralpha * ralpha));
#    ifdef USE_DOUBLE_PRECISION
    const real erfcAlphaR = erfc(ralpha);
#    else
    // This approximation for erfc is from Abramowitz and Stegun (1964) p. 299.
    // They cite the following as the original source: C. Hastings, Jr.,
    // Approximations for Digital Computers (1955).  It has a maximum error
    // of 1.5e-7.

    const real t = RECIP(1.0f + 0.3275911f * ralpha);
    const real erfcAlphaR
      = (0.254829592f + (-0.284496736f + (1.421413741f + (-1.453152027f + 1.061405429f * t) * t) * t) * t) * t * exp2a;
#    endif
    real bn0 = erfcAlphaR * rI;
    real alsq2 = 2 * EWALD_ALPHA * EWALD_ALPHA;
    real alsq2n = RECIP(SQRT_PI * EWALD_ALPHA);
    alsq2n *= alsq2;
    real bn1 = (bn0 + alsq2n * exp2a) * rI2;

    alsq2n *= alsq2;
    real bn2 = (3 * bn1 + alsq2n * exp2a) * rI2;

    // compute the error function scaled and unscaled terms

    real damp = atom1.damp * atom2.damp;
    real ratio = (r / damp);
    ratio = ratio * ratio * ratio;
    float pgamma = atom1.thole < atom2.thole ? atom1.thole : atom2.thole;
    damp = damp == 0 ? 0 : -pgamma * ratio;
    real expdamp = EXP(damp);
    real dsc3 = 1 - expdamp;
    real dsc5 = 1 - expdamp * (1 - damp);
    real r3 = (r * r2);
    real r5 = (r3 * r2);
    real rr3 = (1 - dsc3) / r3;
    real rr5 = 3 * (1 - dsc5) / r5;

    scale1 = rr3 - bn1;
    scale2 = bn2 - rr5;
  } else {
    scale1 = 0;
    scale2 = 0;
  }
  real dDotDelta = scale2 * dot(deltaR, atom2.inducedDipole);
  atom1.field += scale1 * atom2.inducedDipole + dDotDelta * deltaR;
  dDotDelta = scale2 * dot(deltaR, atom2.inducedDipolePolar);
  atom1.fieldPolar += scale1 * atom2.inducedDipolePolar + dDotDelta * deltaR;
  dDotDelta = scale2 * dot(deltaR, atom1.inducedDipole);
  atom2.field += scale1 * atom1.inducedDipole + dDotDelta * deltaR;
  dDotDelta = scale2 * dot(deltaR, atom1.inducedDipolePolar);
  atom2.fieldPolar += scale1 * atom1.inducedDipolePolar + dDotDelta * deltaR;
}

// see also: void computeOneInteraction() in multipoleInducedField.cu
__device__ void tcgPairField1(AtomData1& atom1, AtomData1& atom2, real3 deltaR, bool isSelfInteraction)
{
  if (isSelfInteraction)
    return;
  real scale1, scale2, scale3;
  real r2 = dot(deltaR, deltaR);
  if (r2 < CUTOFF_SQUARED) {
    real rI = RSQRT(r2);
    real r = RECIP(rI);
    real rI2 = rI * rI;

    // calculate the error function damping terms

    real ralpha = EWALD_ALPHA * r;
    real exp2a = EXP(-(ralpha * ralpha));
#    ifdef USE_DOUBLE_PRECISION
    const real erfcAlphaR = erfc(ralpha);
#    else
    // This approximation for erfc is from Abramowitz and Stegun (1964) p. 299.
    // They cite the following as the original source: C. Hastings, Jr.,
    // Approximations for Digital Computers (1955).  It has a maximum error
    // of 1.5e-7.

    const real t = RECIP(1.0f + 0.3275911f * ralpha);
    const real erfcAlphaR
      = (0.254829592f + (-0.284496736f + (1.421413741f + (-1.453152027f + 1.061405429f * t) * t) * t) * t) * t * exp2a;
#    endif
    real bn0 = erfcAlphaR * rI;
    real alsq2 = 2 * EWALD_ALPHA * EWALD_ALPHA;
    real alsq2n = RECIP(SQRT_PI * EWALD_ALPHA);
    alsq2n *= alsq2;
    real bn1 = (bn0 + alsq2n * exp2a) * rI2;

    alsq2n *= alsq2;
    real bn2 = (3 * bn1 + alsq2n * exp2a) * rI2;

    alsq2n *= alsq2;
    real bn3 = (5 * bn2 + alsq2n * exp2a) * rI2;

    // compute the error function scaled and unscaled terms

    real damp = atom1.damp * atom2.damp;
    real ratio = (r / damp);
    ratio = ratio * ratio * ratio;
    float pgamma = atom1.thole < atom2.thole ? atom1.thole : atom2.thole;
    damp = damp == 0 ? 0 : -pgamma * ratio;
    real expdamp = EXP(damp);
    real dsc3 = 1 - expdamp;
    real dsc5 = 1 - expdamp * (1 - damp);
    real dsc7 = 1 - (1 - damp + (0.6f * damp * damp)) * expdamp;
    real r3 = (r * r2);
    real r5 = (r3 * r2);
    real r7 = (r5 * r2);
    real rr3 = (1 - dsc3) / r3;
    real rr5 = 3 * (1 - dsc5) / r5;
    real rr7 = 15 * (1 - dsc7) / r7;

    scale1 = rr3 - bn1;
    scale2 = bn2 - rr5;
    scale3 = bn3 - rr7;
  } else {
    scale1 = 0;
    scale2 = 0;
    scale3 = 0;
  }
  real dDotDelta = scale2 * dot(deltaR, atom2.inducedDipole);
  atom1.field += scale1 * atom2.inducedDipole + dDotDelta * deltaR;
  dDotDelta = scale2 * dot(deltaR, atom2.inducedDipolePolar);
  atom1.fieldPolar += scale1 * atom2.inducedDipolePolar + dDotDelta * deltaR;
  dDotDelta = scale2 * dot(deltaR, atom1.inducedDipole);
  atom2.field += scale1 * atom1.inducedDipole + dDotDelta * deltaR;
  dDotDelta = scale2 * dot(deltaR, atom1.inducedDipolePolar);
  atom2.fieldPolar += scale1 * atom1.inducedDipolePolar + dDotDelta * deltaR;
  // IF CALC GRADIENT
  // Compute and store the field gradients for later use.

  real3 dipole = atom1.inducedDipole;
  real muDotR = dipole.x * deltaR.x + dipole.y * deltaR.y + dipole.z * deltaR.z;
  atom2.fieldGradient[0] -= muDotR * deltaR.x * deltaR.x * scale3 - (2 * dipole.x * deltaR.x + muDotR) * scale2;
  atom2.fieldGradient[1] -= muDotR * deltaR.y * deltaR.y * scale3 - (2 * dipole.y * deltaR.y + muDotR) * scale2;
  atom2.fieldGradient[2] -= muDotR * deltaR.z * deltaR.z * scale3 - (2 * dipole.z * deltaR.z + muDotR) * scale2;
  atom2.fieldGradient[3]
    -= muDotR * deltaR.x * deltaR.y * scale3 - (dipole.x * deltaR.y + dipole.y * deltaR.x) * scale2;
  atom2.fieldGradient[4]
    -= muDotR * deltaR.x * deltaR.z * scale3 - (dipole.x * deltaR.z + dipole.z * deltaR.x) * scale2;
  atom2.fieldGradient[5]
    -= muDotR * deltaR.y * deltaR.z * scale3 - (dipole.y * deltaR.z + dipole.z * deltaR.y) * scale2;

  dipole = atom1.inducedDipolePolar;
  muDotR = dipole.x * deltaR.x + dipole.y * deltaR.y + dipole.z * deltaR.z;
  atom2.fieldGradientPolar[0] -= muDotR * deltaR.x * deltaR.x * scale3 - (2 * dipole.x * deltaR.x + muDotR) * scale2;
  atom2.fieldGradientPolar[1] -= muDotR * deltaR.y * deltaR.y * scale3 - (2 * dipole.y * deltaR.y + muDotR) * scale2;
  atom2.fieldGradientPolar[2] -= muDotR * deltaR.z * deltaR.z * scale3 - (2 * dipole.z * deltaR.z + muDotR) * scale2;
  atom2.fieldGradientPolar[3]
    -= muDotR * deltaR.x * deltaR.y * scale3 - (dipole.x * deltaR.y + dipole.y * deltaR.x) * scale2;
  atom2.fieldGradientPolar[4]
    -= muDotR * deltaR.x * deltaR.z * scale3 - (dipole.x * deltaR.z + dipole.z * deltaR.x) * scale2;
  atom2.fieldGradientPolar[5]
    -= muDotR * deltaR.y * deltaR.z * scale3 - (dipole.y * deltaR.z + dipole.z * deltaR.y) * scale2;

  dipole = atom2.inducedDipole;
  muDotR = dipole.x * deltaR.x + dipole.y * deltaR.y + dipole.z * deltaR.z;
  atom1.fieldGradient[0] += muDotR * deltaR.x * deltaR.x * scale3 - (2 * dipole.x * deltaR.x + muDotR) * scale2;
  atom1.fieldGradient[1] += muDotR * deltaR.y * deltaR.y * scale3 - (2 * dipole.y * deltaR.y + muDotR) * scale2;
  atom1.fieldGradient[2] += muDotR * deltaR.z * deltaR.z * scale3 - (2 * dipole.z * deltaR.z + muDotR) * scale2;
  atom1.fieldGradient[3]
    += muDotR * deltaR.x * deltaR.y * scale3 - (dipole.x * deltaR.y + dipole.y * deltaR.x) * scale2;
  atom1.fieldGradient[4]
    += muDotR * deltaR.x * deltaR.z * scale3 - (dipole.x * deltaR.z + dipole.z * deltaR.x) * scale2;
  atom1.fieldGradient[5]
    += muDotR * deltaR.y * deltaR.z * scale3 - (dipole.y * deltaR.z + dipole.z * deltaR.y) * scale2;

  dipole = atom2.inducedDipolePolar;
  muDotR = dipole.x * deltaR.x + dipole.y * deltaR.y + dipole.z * deltaR.z;
  atom1.fieldGradientPolar[0] += muDotR * deltaR.x * deltaR.x * scale3 - (2 * dipole.x * deltaR.x + muDotR) * scale2;
  atom1.fieldGradientPolar[1] += muDotR * deltaR.y * deltaR.y * scale3 - (2 * dipole.y * deltaR.y + muDotR) * scale2;
  atom1.fieldGradientPolar[2] += muDotR * deltaR.z * deltaR.z * scale3 - (2 * dipole.z * deltaR.z + muDotR) * scale2;
  atom1.fieldGradientPolar[3]
    += muDotR * deltaR.x * deltaR.y * scale3 - (dipole.x * deltaR.y + dipole.y * deltaR.x) * scale2;
  atom1.fieldGradientPolar[4]
    += muDotR * deltaR.x * deltaR.z * scale3 - (dipole.x * deltaR.z + dipole.z * deltaR.x) * scale2;
  atom1.fieldGradientPolar[5]
    += muDotR * deltaR.y * deltaR.z * scale3 - (dipole.y * deltaR.z + dipole.z * deltaR.y) * scale2;
  // END IF
}
#  else  // USE_EWALD
#  endif // USE_EWALD

// see also: void computeInducedField(...) in multipoleInducedField.cu
extern "C" __global__ void tcgUField0(unsigned long long* __restrict__ field,
  unsigned long long* __restrict__ fieldPolar, const real4* __restrict__ posq,
  const ushort2* __restrict__ exclusionTiles, const real* __restrict__ inducedDipole,
  const real* __restrict__ inducedDipolePolar, unsigned int startTileIndex, unsigned int numTileIndices,
  // IF CALC GRADIENT
  // unsigned long long* __restrict__ fieldGradient, unsigned long long* __restrict__ fieldGradientPolar,
  // END IF
  // #ifdef USE_CUTOFF
  const int* __restrict__ tiles, const unsigned int* __restrict__ interactionCount, real4 periodicBoxSize,
  real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles,
  const real4* __restrict__ blockCenter, const unsigned int* __restrict__ interactingAtoms,
  // #endif
  const float2* __restrict__ dampingAndThole)
{
  const unsigned int totalWarps = (blockDim.x * gridDim.x) / TILE_SIZE;
  const unsigned int warp = (blockIdx.x * blockDim.x + threadIdx.x) / TILE_SIZE;
  const unsigned int tgx = threadIdx.x & (TILE_SIZE - 1);
  const unsigned int tbx = threadIdx.x - tgx;
  __shared__ AtomData0 localData[THREAD_BLOCK_SIZE_0];
  __shared__ int atomIndices[THREAD_BLOCK_SIZE_0];
  __shared__ volatile int skipTiles[THREAD_BLOCK_SIZE_0];

  // First loop: process tiles that contain exclusions.

  const unsigned int firstExclusionTile
    = FIRST_EXCLUSION_TILE + warp * (LAST_EXCLUSION_TILE - FIRST_EXCLUSION_TILE) / totalWarps;
  const unsigned int lastExclusionTile
    = FIRST_EXCLUSION_TILE + (warp + 1) * (LAST_EXCLUSION_TILE - FIRST_EXCLUSION_TILE) / totalWarps;
  for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
    const ushort2 tileIndices = exclusionTiles[pos];
    const unsigned int x = tileIndices.x;
    const unsigned int y = tileIndices.y;
    AtomData0 data;
    zeroAtomData0(data);
    unsigned int atom1 = x * TILE_SIZE + tgx;
    loadAtomData0(data, atom1, posq, inducedDipole, inducedDipolePolar, dampingAndThole);

    if (x == y) {
      // This tile is on the diagonal.
      localData[threadIdx.x].pos = data.pos;
      localData[threadIdx.x].inducedDipole = data.inducedDipole;
      localData[threadIdx.x].inducedDipolePolar = data.inducedDipolePolar;
      localData[threadIdx.x].thole = data.thole;
      localData[threadIdx.x].damp = data.damp;
      for (unsigned int j = 0; j < TILE_SIZE; j++) {
        real3 delta = localData[tbx + j].pos - data.pos;
#  ifdef USE_PERIODIC
        APPLY_PERIODIC_TO_DELTA(delta)
#  endif
        int atom2 = y * TILE_SIZE + j;
        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)
          tcgPairField0(data, localData[tbx + j], delta, atom1 == atom2);
      }
    } else {
      // This is an off-diagonal tile.
      loadAtomData0(
        localData[threadIdx.x], y * TILE_SIZE + tgx, posq, inducedDipole, inducedDipolePolar, dampingAndThole);
      zeroAtomData0(localData[threadIdx.x]);
      unsigned int tj = tgx;
      for (unsigned int j = 0; j < TILE_SIZE; j++) {
        real3 delta = localData[tbx + tj].pos - data.pos;
#  ifdef USE_PERIODIC
        APPLY_PERIODIC_TO_DELTA(delta)
#  endif
        int atom2 = y * TILE_SIZE + j;
        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)
          tcgPairField0(data, localData[tbx + tj], delta, false);
        tj = (tj + 1) & (TILE_SIZE - 1);
      }
    }

    // Write results.

    unsigned int offset = x * TILE_SIZE + tgx;
    SAVE_ATOM_DATA0(offset, data)
    if (x != y) {
      offset = y * TILE_SIZE + tgx;
      SAVE_ATOM_DATA0(offset, localData[threadIdx.x])
    }
  }

  // Second loop: tiles without exclusions, either from the neighbor list (with
  // cutoff) or just enumerating all of them (no cutoff).

#  ifdef USE_CUTOFF
  const unsigned int numTiles = interactionCount[0];
  if (numTiles > maxTiles)
    return; // There wasn't enough memory for the neighbor list.
  int pos = (int)(numTiles > maxTiles ? startTileIndex + warp * (long long)numTileIndices / totalWarps
                                      : warp * (long long)numTiles / totalWarps);
  int end = (int)(numTiles > maxTiles ? startTileIndex + (warp + 1) * (long long)numTileIndices / totalWarps
                                      : (warp + 1) * (long long)numTiles / totalWarps);
#  else
  const unsigned int numTiles = numTileIndices;
  int pos = (int)(startTileIndex + warp * (long long)numTiles / totalWarps);
  int end = (int)(startTileIndex + (warp + 1) * (long long)numTiles / totalWarps);
#  endif
  int skipBase = 0;
  int currentSkipIndex = tbx;
  skipTiles[threadIdx.x] = -1;

  while (pos < end) {
    bool includeTile = true;

    // Extract the coordinates of this tile.

    int x, y;
#  ifdef USE_CUTOFF
    x = tiles[pos];
#  else
    y = (int)floor(NUM_BLOCKS + 0.5f - SQRT((NUM_BLOCKS + 0.5f) * (NUM_BLOCKS + 0.5f) - 2 * pos));
    x = (pos - y * NUM_BLOCKS + y * (y + 1) / 2);
    if (x < y || x >= NUM_BLOCKS) {
      // Occasionally happens due to roundoff error.
      y += (x < y ? -1 : 1);
      x = (pos - y * NUM_BLOCKS + y * (y + 1) / 2);
    }

    // Skip over tiles that have exclusions, since they were already
    // processed.

    while (skipTiles[tbx + TILE_SIZE - 1] < pos) {
      if (skipBase + tgx < NUM_TILES_WITH_EXCLUSIONS) {
        ushort2 tile = exclusionTiles[skipBase + tgx];
        skipTiles[threadIdx.x] = tile.x + tile.y * NUM_BLOCKS - tile.y * (tile.y + 1) / 2;
      } else
        skipTiles[threadIdx.x] = end;
      skipBase += TILE_SIZE;
      currentSkipIndex = tbx;
    }
    while (skipTiles[currentSkipIndex] < pos)
      currentSkipIndex++;
    includeTile = (skipTiles[currentSkipIndex] != pos);
#  endif
    if (includeTile) {
      unsigned int atom1 = x * TILE_SIZE + tgx;

      // Load atom data for this tile.

      AtomData0 data;
      zeroAtomData0(data);
      loadAtomData0(data, atom1, posq, inducedDipole, inducedDipolePolar, dampingAndThole);
#  ifdef USE_CUTOFF
      unsigned int j = interactingAtoms[pos * TILE_SIZE + tgx];
#  else
      unsigned int j = y * TILE_SIZE + tgx;
#  endif
      atomIndices[threadIdx.x] = j;
      loadAtomData0(localData[threadIdx.x], j, posq, inducedDipole, inducedDipolePolar, dampingAndThole);
      zeroAtomData0(localData[threadIdx.x]);

      // Compute the full set of interactions in this tile.

      unsigned int tj = tgx;
      for (j = 0; j < TILE_SIZE; j++) {
        real3 delta = localData[tbx + tj].pos - data.pos;
#  ifdef USE_PERIODIC
        APPLY_PERIODIC_TO_DELTA(delta)
#  endif
        int atom2 = atomIndices[tbx + tj];
        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)
          tcgPairField0(data, localData[tbx + tj], delta, false);
        tj = (tj + 1) & (TILE_SIZE - 1);
      }

      // Write results.

      unsigned int offset = x * TILE_SIZE + tgx;
      SAVE_ATOM_DATA0(offset, data)
#  ifdef USE_CUTOFF
      offset = atomIndices[threadIdx.x];
#  else
      offset = y * TILE_SIZE + tgx;
#  endif
      SAVE_ATOM_DATA0(offset, localData[threadIdx.x])
    }
    pos++;
  }
}

// see also: void computeInducedField(...) in multipoleInducedField.cu
extern "C" __global__ void tcgUField1(unsigned long long* __restrict__ field,
  unsigned long long* __restrict__ fieldPolar, const real4* __restrict__ posq,
  const ushort2* __restrict__ exclusionTiles, const real* __restrict__ inducedDipole,
  const real* __restrict__ inducedDipolePolar, unsigned int startTileIndex, unsigned int numTileIndices,
  // IF CALC GRADIENT
  unsigned long long* __restrict__ fieldGradient, unsigned long long* __restrict__ fieldGradientPolar,
  // END IF
  // #ifdef USE_CUTOFF
  const int* __restrict__ tiles, const unsigned int* __restrict__ interactionCount, real4 periodicBoxSize,
  real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles,
  const real4* __restrict__ blockCenter, const unsigned int* __restrict__ interactingAtoms,
  // #endif
  const float2* __restrict__ dampingAndThole)
{
  const unsigned int totalWarps = (blockDim.x * gridDim.x) / TILE_SIZE;
  const unsigned int warp = (blockIdx.x * blockDim.x + threadIdx.x) / TILE_SIZE;
  const unsigned int tgx = threadIdx.x & (TILE_SIZE - 1);
  const unsigned int tbx = threadIdx.x - tgx;
  __shared__ AtomData1 localData[THREAD_BLOCK_SIZE_1];
  __shared__ int atomIndices[THREAD_BLOCK_SIZE_1];
  __shared__ volatile int skipTiles[THREAD_BLOCK_SIZE_1];

  // First loop: process tiles that contain exclusions.

  const unsigned int firstExclusionTile
    = FIRST_EXCLUSION_TILE + warp * (LAST_EXCLUSION_TILE - FIRST_EXCLUSION_TILE) / totalWarps;
  const unsigned int lastExclusionTile
    = FIRST_EXCLUSION_TILE + (warp + 1) * (LAST_EXCLUSION_TILE - FIRST_EXCLUSION_TILE) / totalWarps;
  for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
    const ushort2 tileIndices = exclusionTiles[pos];
    const unsigned int x = tileIndices.x;
    const unsigned int y = tileIndices.y;
    AtomData1 data;
    zeroAtomData1(data);
    unsigned int atom1 = x * TILE_SIZE + tgx;
    loadAtomData1(data, atom1, posq, inducedDipole, inducedDipolePolar, dampingAndThole);

    if (x == y) {
      // This tile is on the diagonal.
      localData[threadIdx.x].pos = data.pos;
      localData[threadIdx.x].inducedDipole = data.inducedDipole;
      localData[threadIdx.x].inducedDipolePolar = data.inducedDipolePolar;
      localData[threadIdx.x].thole = data.thole;
      localData[threadIdx.x].damp = data.damp;
      for (unsigned int j = 0; j < TILE_SIZE; j++) {
        real3 delta = localData[tbx + j].pos - data.pos;
#  ifdef USE_PERIODIC
        APPLY_PERIODIC_TO_DELTA(delta)
#  endif
        int atom2 = y * TILE_SIZE + j;
        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)
          tcgPairField1(data, localData[tbx + j], delta, atom1 == atom2);
      }
    } else {
      // This is an off-diagonal tile.
      loadAtomData1(
        localData[threadIdx.x], y * TILE_SIZE + tgx, posq, inducedDipole, inducedDipolePolar, dampingAndThole);
      zeroAtomData1(localData[threadIdx.x]);
      unsigned int tj = tgx;
      for (unsigned int j = 0; j < TILE_SIZE; j++) {
        real3 delta = localData[tbx + tj].pos - data.pos;
#  ifdef USE_PERIODIC
        APPLY_PERIODIC_TO_DELTA(delta)
#  endif
        int atom2 = y * TILE_SIZE + j;
        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)
          tcgPairField1(data, localData[tbx + tj], delta, false);
        tj = (tj + 1) & (TILE_SIZE - 1);
      }
    }

    // Write results.

    unsigned int offset = x * TILE_SIZE + tgx;
    SAVE_ATOM_DATA1(offset, data)
    if (x != y) {
      offset = y * TILE_SIZE + tgx;
      SAVE_ATOM_DATA1(offset, localData[threadIdx.x])
    }
  }

  // Second loop: tiles without exclusions, either from the neighbor list
  // (with cutoff) or just enumerating all of them (no cutoff).

#  ifdef USE_CUTOFF
  const unsigned int numTiles = interactionCount[0];
  if (numTiles > maxTiles)
    return; // There wasn't enough memory for the neighbor list.
  int pos = (int)(numTiles > maxTiles ? startTileIndex + warp * (long long)numTileIndices / totalWarps
                                      : warp * (long long)numTiles / totalWarps);
  int end = (int)(numTiles > maxTiles ? startTileIndex + (warp + 1) * (long long)numTileIndices / totalWarps
                                      : (warp + 1) * (long long)numTiles / totalWarps);
#  else
  const unsigned int numTiles = numTileIndices;
  int pos = (int)(startTileIndex + warp * (long long)numTiles / totalWarps);
  int end = (int)(startTileIndex + (warp + 1) * (long long)numTiles / totalWarps);
#  endif
  int skipBase = 0;
  int currentSkipIndex = tbx;
  skipTiles[threadIdx.x] = -1;

  while (pos < end) {
    bool includeTile = true;

    // Extract the coordinates of this tile.

    int x, y;
#  ifdef USE_CUTOFF
    x = tiles[pos];
#  else
    y = (int)floor(NUM_BLOCKS + 0.5f - SQRT((NUM_BLOCKS + 0.5f) * (NUM_BLOCKS + 0.5f) - 2 * pos));
    x = (pos - y * NUM_BLOCKS + y * (y + 1) / 2);
    if (x < y || x >= NUM_BLOCKS) {
      // Occasionally happens due to roundoff error.
      y += (x < y ? -1 : 1);
      x = (pos - y * NUM_BLOCKS + y * (y + 1) / 2);
    }

    // Skip over tiles that have exclusions, since they were already
    // processed.

    while (skipTiles[tbx + TILE_SIZE - 1] < pos) {
      if (skipBase + tgx < NUM_TILES_WITH_EXCLUSIONS) {
        ushort2 tile = exclusionTiles[skipBase + tgx];
        skipTiles[threadIdx.x] = tile.x + tile.y * NUM_BLOCKS - tile.y * (tile.y + 1) / 2;
      } else
        skipTiles[threadIdx.x] = end;
      skipBase += TILE_SIZE;
      currentSkipIndex = tbx;
    }
    while (skipTiles[currentSkipIndex] < pos)
      currentSkipIndex++;
    includeTile = (skipTiles[currentSkipIndex] != pos);
#  endif
    if (includeTile) {
      unsigned int atom1 = x * TILE_SIZE + tgx;

      // Load atom data for this tile.

      AtomData1 data;
      zeroAtomData1(data);
      loadAtomData1(data, atom1, posq, inducedDipole, inducedDipolePolar, dampingAndThole);
#  ifdef USE_CUTOFF
      unsigned int j = interactingAtoms[pos * TILE_SIZE + tgx];
#  else
      unsigned int j = y * TILE_SIZE + tgx;
#  endif
      atomIndices[threadIdx.x] = j;
      loadAtomData1(localData[threadIdx.x], j, posq, inducedDipole, inducedDipolePolar, dampingAndThole);
      zeroAtomData1(localData[threadIdx.x]);

      // Compute the full set of interactions in this tile.

      unsigned int tj = tgx;
      for (j = 0; j < TILE_SIZE; j++) {
        real3 delta = localData[tbx + tj].pos - data.pos;
#  ifdef USE_PERIODIC
        APPLY_PERIODIC_TO_DELTA(delta)
#  endif
        int atom2 = atomIndices[tbx + tj];
        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)
          tcgPairField1(data, localData[tbx + tj], delta, false);
        tj = (tj + 1) & (TILE_SIZE - 1);
      }

      // Write results.

      unsigned int offset = x * TILE_SIZE + tgx;
      SAVE_ATOM_DATA1(offset, data)
#  ifdef USE_CUTOFF
      offset = atomIndices[threadIdx.x];
#  else
      offset = y * TILE_SIZE + tgx;
#  endif
      SAVE_ATOM_DATA1(offset, localData[threadIdx.x])
    }
    pos++;
  }
}

////////////////////////////////////////////////////////////////////////////////

// see also: void addExtrapolatedFieldGradientToForce() in multipoleInducedField.cu
extern "C" __global__ void tcgAddMutualForce(long long* __restrict__ forceBuffers, real* __restrict__ iuad,
  real* __restrict__ iuap, real* __restrict__ iubd, real* __restrict__ iubp, long long* __restrict__ iuadfg,
  long long* __restrict__ iuapfg, long long* __restrict__ iubdfg, long long* __restrict__ iubpfg)
{
  const real scale = 0.5f * ENERGY_SCALE_FACTOR;
  for (int atom = blockIdx.x * blockDim.x + threadIdx.x; atom < NUM_ATOMS; atom += blockDim.x * gridDim.x) {
    real fx = 0, fy = 0, fz = 0;
    for (int l = 0; l < M_TCGNAB; ++l) {
      int index1 = 3 * (l * PADDED_NUM_ATOMS + atom);
      int index2 = 6 * (l * PADDED_NUM_ATOMS + atom);

      real uad[] = {iuad[index1], iuad[index1 + 1], iuad[index1 + 2]};
      real uap[] = {iuap[index1], iuap[index1 + 1], iuap[index1 + 2]};
      real ubd[] = {iubd[index1], iubd[index1 + 1], iubd[index1 + 2]};
      real ubp[] = {iubp[index1], iubp[index1 + 1], iubp[index1 + 2]};
      long long uadfg[] = {iuadfg[index2], iuadfg[index2 + 1], iuadfg[index2 + 2], iuadfg[index2 + 3],
        iuadfg[index2 + 4], iuadfg[index2 + 5]};
      long long ubdfg[] = {iubdfg[index2], iubdfg[index2 + 1], iubdfg[index2 + 2], iubdfg[index2 + 3],
        iubdfg[index2 + 4], iubdfg[index2 + 5]};
      long long uapfg[] = {iuapfg[index2], iuapfg[index2 + 1], iuapfg[index2 + 2], iuapfg[index2 + 3],
        iuapfg[index2 + 4], iuapfg[index2 + 5]};
      long long ubpfg[] = {iubpfg[index2], iubpfg[index2 + 1], iubpfg[index2 + 2], iubpfg[index2 + 3],
        iubpfg[index2 + 4], iubpfg[index2 + 5]};

      fx += scale
        * (uad[0] * ubpfg[0] + uad[1] * ubpfg[3] + uad[2] * ubpfg[4] + uap[0] * ubdfg[0] + uap[1] * ubdfg[3]
            + uap[2] * ubdfg[4] + ubd[0] * uapfg[0] + ubd[1] * uapfg[3] + ubd[2] * uapfg[4] + ubp[0] * uadfg[0]
            + ubp[1] * uadfg[3] + ubp[2] * uadfg[4]);
      fy += scale
        * (uad[0] * ubpfg[3] + uad[1] * ubpfg[1] + uad[2] * ubpfg[5] + uap[0] * ubdfg[3] + uap[1] * ubdfg[1]
            + uap[2] * ubdfg[5] + ubd[0] * uapfg[3] + ubd[1] * uapfg[1] + ubd[2] * uapfg[5] + ubp[0] * uadfg[3]
            + ubp[1] * uadfg[1] + ubp[2] * uadfg[5]);
      fz += scale
        * (uad[0] * ubpfg[4] + uad[1] * ubpfg[5] + uad[2] * ubpfg[2] + uap[0] * ubdfg[4] + uap[1] * ubdfg[5]
            + uap[2] * ubdfg[2] + ubd[0] * uapfg[4] + ubd[1] * uapfg[5] + ubd[2] * uapfg[2] + ubp[0] * uadfg[4]
            + ubp[1] * uadfg[5] + ubp[2] * uadfg[2]);
    }

    forceBuffers[atom] += (long long)fx;
    forceBuffers[atom + PADDED_NUM_ATOMS] += (long long)fy;
    forceBuffers[atom + PADDED_NUM_ATOMS * 2] += (long long)fz;
  }
}

// see also: void recordInducedFieldDipoles(...) in multipolePme.cu
extern "C" __global__ void tcgRecordPmeInducedField0(const real* __restrict__ phid, real* const __restrict__ phip,
  long long* __restrict__ inducedField, long long* __restrict__ inducedFieldPolar,
  const real* __restrict__ inducedDipole, const real* __restrict__ inducedDipolePolar,
  // IF CALC GRADIENT
  // unsigned long long* __restrict__ fieldGradient, unsigned long long* __restrict__ fieldGradientPolar,
  // END IF
  real3 recipBoxVecX, real3 recipBoxVecY, real3 recipBoxVecZ)
{
  __shared__ real fracToCart[3][3];
  if (threadIdx.x == 0) {
    fracToCart[0][0] = GRID_SIZE_X * recipBoxVecX.x;
    fracToCart[1][0] = GRID_SIZE_X * recipBoxVecY.x;
    fracToCart[2][0] = GRID_SIZE_X * recipBoxVecZ.x;
    fracToCart[0][1] = GRID_SIZE_Y * recipBoxVecX.y;
    fracToCart[1][1] = GRID_SIZE_Y * recipBoxVecY.y;
    fracToCart[2][1] = GRID_SIZE_Y * recipBoxVecZ.y;
    fracToCart[0][2] = GRID_SIZE_Z * recipBoxVecX.z;
    fracToCart[1][2] = GRID_SIZE_Z * recipBoxVecY.z;
    fracToCart[2][2] = GRID_SIZE_Z * recipBoxVecZ.z;
  }
  __syncthreads();
  const real selfDipoleScale = (4 / (real)3) * (EWALD_ALPHA * EWALD_ALPHA * EWALD_ALPHA) / SQRT_PI;
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < NUM_ATOMS; i += blockDim.x * gridDim.x) {
    inducedField[i] -= (long long)(0x100000000
      * (phid[i + NUM_ATOMS] * fracToCart[0][0] + phid[i + NUM_ATOMS * 2] * fracToCart[0][1]
          + phid[i + NUM_ATOMS * 3] * fracToCart[0][2] - selfDipoleScale * inducedDipole[3 * i]));
    inducedField[i + PADDED_NUM_ATOMS] -= (long long)(0x100000000
      * (phid[i + NUM_ATOMS] * fracToCart[1][0] + phid[i + NUM_ATOMS * 2] * fracToCart[1][1]
          + phid[i + NUM_ATOMS * 3] * fracToCart[1][2] - selfDipoleScale * inducedDipole[3 * i + 1]));
    inducedField[i + PADDED_NUM_ATOMS * 2] -= (long long)(0x100000000
      * (phid[i + NUM_ATOMS] * fracToCart[2][0] + phid[i + NUM_ATOMS * 2] * fracToCart[2][1]
          + phid[i + NUM_ATOMS * 3] * fracToCart[2][2] - selfDipoleScale * inducedDipole[3 * i + 2]));
    inducedFieldPolar[i] -= (long long)(0x100000000
      * (phip[i + NUM_ATOMS] * fracToCart[0][0] + phip[i + NUM_ATOMS * 2] * fracToCart[0][1]
          + phip[i + NUM_ATOMS * 3] * fracToCart[0][2] - selfDipoleScale * inducedDipolePolar[3 * i]));
    inducedFieldPolar[i + PADDED_NUM_ATOMS] -= (long long)(0x100000000
      * (phip[i + NUM_ATOMS] * fracToCart[1][0] + phip[i + NUM_ATOMS * 2] * fracToCart[1][1]
          + phip[i + NUM_ATOMS * 3] * fracToCart[1][2] - selfDipoleScale * inducedDipolePolar[3 * i + 1]));
    inducedFieldPolar[i + PADDED_NUM_ATOMS * 2] -= (long long)(0x100000000
      * (phip[i + NUM_ATOMS] * fracToCart[2][0] + phip[i + NUM_ATOMS * 2] * fracToCart[2][1]
          + phip[i + NUM_ATOMS * 3] * fracToCart[2][2] - selfDipoleScale * inducedDipolePolar[3 * i + 2]));
  }
}

// see also: void recordInducedFieldDipoles(...) in multipolePme.cu
extern "C" __global__ void tcgRecordPmeInducedField1(const real* __restrict__ phid, real* const __restrict__ phip,
  long long* __restrict__ inducedField, long long* __restrict__ inducedFieldPolar,
  const real* __restrict__ inducedDipole, const real* __restrict__ inducedDipolePolar,
  // IF CALC GRADIENT
  unsigned long long* __restrict__ fieldGradient, unsigned long long* __restrict__ fieldGradientPolar,
  // END IF
  real3 recipBoxVecX, real3 recipBoxVecY, real3 recipBoxVecZ)
{
  __shared__ real fracToCart[3][3];
  if (threadIdx.x == 0) {
    fracToCart[0][0] = GRID_SIZE_X * recipBoxVecX.x;
    fracToCart[1][0] = GRID_SIZE_X * recipBoxVecY.x;
    fracToCart[2][0] = GRID_SIZE_X * recipBoxVecZ.x;
    fracToCart[0][1] = GRID_SIZE_Y * recipBoxVecX.y;
    fracToCart[1][1] = GRID_SIZE_Y * recipBoxVecY.y;
    fracToCart[2][1] = GRID_SIZE_Y * recipBoxVecZ.y;
    fracToCart[0][2] = GRID_SIZE_Z * recipBoxVecX.z;
    fracToCart[1][2] = GRID_SIZE_Z * recipBoxVecY.z;
    fracToCart[2][2] = GRID_SIZE_Z * recipBoxVecZ.z;
  }
  __syncthreads();

  const real selfDipoleScale = (4 / (real)3) * (EWALD_ALPHA * EWALD_ALPHA * EWALD_ALPHA) / SQRT_PI;
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < NUM_ATOMS; i += blockDim.x * gridDim.x) {
    inducedField[i] -= (long long)(0x100000000
      * (phid[i + NUM_ATOMS] * fracToCart[0][0] + phid[i + NUM_ATOMS * 2] * fracToCart[0][1]
          + phid[i + NUM_ATOMS * 3] * fracToCart[0][2] - selfDipoleScale * inducedDipole[3 * i]));
    inducedField[i + PADDED_NUM_ATOMS] -= (long long)(0x100000000
      * (phid[i + NUM_ATOMS] * fracToCart[1][0] + phid[i + NUM_ATOMS * 2] * fracToCart[1][1]
          + phid[i + NUM_ATOMS * 3] * fracToCart[1][2] - selfDipoleScale * inducedDipole[3 * i + 1]));
    inducedField[i + PADDED_NUM_ATOMS * 2] -= (long long)(0x100000000
      * (phid[i + NUM_ATOMS] * fracToCart[2][0] + phid[i + NUM_ATOMS * 2] * fracToCart[2][1]
          + phid[i + NUM_ATOMS * 3] * fracToCart[2][2] - selfDipoleScale * inducedDipole[3 * i + 2]));
    inducedFieldPolar[i] -= (long long)(0x100000000
      * (phip[i + NUM_ATOMS] * fracToCart[0][0] + phip[i + NUM_ATOMS * 2] * fracToCart[0][1]
          + phip[i + NUM_ATOMS * 3] * fracToCart[0][2] - selfDipoleScale * inducedDipolePolar[3 * i]));
    inducedFieldPolar[i + PADDED_NUM_ATOMS] -= (long long)(0x100000000
      * (phip[i + NUM_ATOMS] * fracToCart[1][0] + phip[i + NUM_ATOMS * 2] * fracToCart[1][1]
          + phip[i + NUM_ATOMS * 3] * fracToCart[1][2] - selfDipoleScale * inducedDipolePolar[3 * i + 1]));
    inducedFieldPolar[i + PADDED_NUM_ATOMS * 2] -= (long long)(0x100000000
      * (phip[i + NUM_ATOMS] * fracToCart[2][0] + phip[i + NUM_ATOMS * 2] * fracToCart[2][1]
          + phip[i + NUM_ATOMS * 3] * fracToCart[2][2] - selfDipoleScale * inducedDipolePolar[3 * i + 2]));

    // Compute and store the field gradients for later use.

    real EmatD[3][3] = {{phid[i + NUM_ATOMS * 4], phid[i + NUM_ATOMS * 7], phid[i + NUM_ATOMS * 8]},
      {phid[i + NUM_ATOMS * 7], phid[i + NUM_ATOMS * 5], phid[i + NUM_ATOMS * 9]},
      {phid[i + NUM_ATOMS * 8], phid[i + NUM_ATOMS * 9], phid[i + NUM_ATOMS * 6]}};
    real Exx = 0, Eyy = 0, Ezz = 0, Exy = 0, Exz = 0, Eyz = 0;
    for (int k = 0; k < 3; ++k) {
      for (int l = 0; l < 3; ++l) {
        Exx += fracToCart[0][k] * EmatD[k][l] * fracToCart[0][l];
        Eyy += fracToCart[1][k] * EmatD[k][l] * fracToCart[1][l];
        Ezz += fracToCart[2][k] * EmatD[k][l] * fracToCart[2][l];
        Exy += fracToCart[0][k] * EmatD[k][l] * fracToCart[1][l];
        Exz += fracToCart[0][k] * EmatD[k][l] * fracToCart[2][l];
        Eyz += fracToCart[1][k] * EmatD[k][l] * fracToCart[2][l];
      }
    }
    /*
    atomicAdd(&fieldGradient[6 * i + 0], static_cast<unsigned long long>((long long)(-Exx * 0x100000000)));
    atomicAdd(&fieldGradient[6 * i + 1], static_cast<unsigned long long>((long long)(-Eyy * 0x100000000)));
    atomicAdd(&fieldGradient[6 * i + 2], static_cast<unsigned long long>((long long)(-Ezz * 0x100000000)));
    atomicAdd(&fieldGradient[6 * i + 3], static_cast<unsigned long long>((long long)(-Exy * 0x100000000)));
    atomicAdd(&fieldGradient[6 * i + 4], static_cast<unsigned long long>((long long)(-Exz * 0x100000000)));
    atomicAdd(&fieldGradient[6 * i + 5], static_cast<unsigned long long>((long long)(-Eyz * 0x100000000)));
    */
    fieldGradient[6 * i + 0] += static_cast<unsigned long long>((long long)(-Exx * 0x100000000));
    fieldGradient[6 * i + 1] += static_cast<unsigned long long>((long long)(-Eyy * 0x100000000));
    fieldGradient[6 * i + 2] += static_cast<unsigned long long>((long long)(-Ezz * 0x100000000));
    fieldGradient[6 * i + 3] += static_cast<unsigned long long>((long long)(-Exy * 0x100000000));
    fieldGradient[6 * i + 4] += static_cast<unsigned long long>((long long)(-Exz * 0x100000000));
    fieldGradient[6 * i + 5] += static_cast<unsigned long long>((long long)(-Eyz * 0x100000000));

    real EmatP[3][3] = {{phip[i + NUM_ATOMS * 4], phip[i + NUM_ATOMS * 7], phip[i + NUM_ATOMS * 8]},
      {phip[i + NUM_ATOMS * 7], phip[i + NUM_ATOMS * 5], phip[i + NUM_ATOMS * 9]},
      {phip[i + NUM_ATOMS * 8], phip[i + NUM_ATOMS * 9], phip[i + NUM_ATOMS * 6]}};
    Exx = 0;
    Eyy = 0;
    Ezz = 0;
    Exy = 0;
    Exz = 0;
    Eyz = 0;
    for (int k = 0; k < 3; ++k) {
      for (int l = 0; l < 3; ++l) {
        Exx += fracToCart[0][k] * EmatP[k][l] * fracToCart[0][l];
        Eyy += fracToCart[1][k] * EmatP[k][l] * fracToCart[1][l];
        Ezz += fracToCart[2][k] * EmatP[k][l] * fracToCart[2][l];
        Exy += fracToCart[0][k] * EmatP[k][l] * fracToCart[1][l];
        Exz += fracToCart[0][k] * EmatP[k][l] * fracToCart[2][l];
        Eyz += fracToCart[1][k] * EmatP[k][l] * fracToCart[2][l];
      }
    }
    /*
    atomicAdd(&fieldGradientPolar[6 * i + 0], static_cast<unsigned long long>((long long)(-Exx * 0x100000000)));
    atomicAdd(&fieldGradientPolar[6 * i + 1], static_cast<unsigned long long>((long long)(-Eyy * 0x100000000)));
    atomicAdd(&fieldGradientPolar[6 * i + 2], static_cast<unsigned long long>((long long)(-Ezz * 0x100000000)));
    atomicAdd(&fieldGradientPolar[6 * i + 3], static_cast<unsigned long long>((long long)(-Exy * 0x100000000)));
    atomicAdd(&fieldGradientPolar[6 * i + 4], static_cast<unsigned long long>((long long)(-Exz * 0x100000000)));
    atomicAdd(&fieldGradientPolar[6 * i + 5], static_cast<unsigned long long>((long long)(-Eyz * 0x100000000)));
    */
    fieldGradientPolar[6 * i + 0] += static_cast<unsigned long long>((long long)(-Exx * 0x100000000));
    fieldGradientPolar[6 * i + 1] += static_cast<unsigned long long>((long long)(-Eyy * 0x100000000));
    fieldGradientPolar[6 * i + 2] += static_cast<unsigned long long>((long long)(-Ezz * 0x100000000));
    fieldGradientPolar[6 * i + 3] += static_cast<unsigned long long>((long long)(-Exy * 0x100000000));
    fieldGradientPolar[6 * i + 4] += static_cast<unsigned long long>((long long)(-Exz * 0x100000000));
    fieldGradientPolar[6 * i + 5] += static_cast<unsigned long long>((long long)(-Eyz * 0x100000000));
  }
}

#endif

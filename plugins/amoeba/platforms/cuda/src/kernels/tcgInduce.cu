#ifdef TCG_POLARIZATION

#  if !defined(TCG_USE_DOUBLE_PRECISION) || !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#    define TCG_ATOMIC_ADD_REAL atomicAdd
#  else
#    define TCG_ATOMIC_ADD_REAL tcgAtomicAddDouble
__device__ double tcgAtomicAddDouble(double* address, double val)
{
  unsigned long long int* address_as_ull = (unsigned long long*)address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}
#  endif

/**
 * Copy from (long long) field/p arrays to (real) uind/p arrays, and reshape the data layout.
 *
 * Parameters:
 * long long, dimension(PADDED_NUM_ATOMS,3) :: field,fieldp
 * real, dimension(3,*) :: uind,uinp
 */
extern "C" __global__ void tcgFieldLongLongToReal(const long long* __restrict__ field,
  const long long* __restrict__ fieldp, real* __restrict__ uind, real* __restrict__ uinp)
{
  const real fieldScale = 1 / (real)0x100000000;
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < NUM_ATOMS; i += blockDim.x * gridDim.x) {
    int index = 3 * i;
    uind[index + 0] = fieldScale * field[i];
    uind[index + 1] = fieldScale * field[i + PADDED_NUM_ATOMS];
    uind[index + 2] = fieldScale * field[i + PADDED_NUM_ATOMS * 2];
    uinp[index + 0] = fieldScale * fieldp[i];
    uinp[index + 1] = fieldScale * fieldp[i + PADDED_NUM_ATOMS];
    uinp[index + 2] = fieldScale * fieldp[i + PADDED_NUM_ATOMS * 2];
  }
}

/**
 * This is a helper kernel for T-dot-Vector [output = T.x], where T = 1/alpha + Tu.
 *
 * [-Tu.x = field] is pre-calculated by the induced dipole field kernels,
 * and this kernel will perform [output = (1/alpha).x - field].
 *
 * Parameters:
 * real, dimension(3,*) :: outd,outp [output]
 * float, dimension(*) :: polarity [alpha]
 * real, dimension(3,*) :: uind,uinp [x]
 * long long, dimension(PADDED_NUM_ATOMS,3) :: field,fieldp [-Tu.x]
 */
extern "C" __global__ void tcgInvAlphaPlusTu(real* __restrict__ outd, real* __restrict__ outp,
  const float* __restrict__ polarity, const real* __restrict__ uind, const real* __restrict__ uinp,
  const long long* __restrict__ field, const long long* __restrict__ fieldp)
{
  const real fieldScale = 1 / (real)0x100000000;
  const float polmin = 1.0e-9;
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < NUM_ATOMS; i += blockDim.x * gridDim.x) {
    int index1 = 3 * i;
    const float poli = max(polmin, polarity[i]);
    const real invPoli = RECIP(poli);
    outd[index1 + 0] = uind[index1 + 0] * invPoli - fieldScale * field[i];
    outd[index1 + 1] = uind[index1 + 1] * invPoli - fieldScale * field[i + PADDED_NUM_ATOMS];
    outd[index1 + 2] = uind[index1 + 2] * invPoli - fieldScale * field[i + PADDED_NUM_ATOMS * 2];
    outp[index1 + 0] = uinp[index1 + 0] * invPoli - fieldScale * fieldp[i];
    outp[index1 + 1] = uinp[index1 + 1] * invPoli - fieldScale * fieldp[i + PADDED_NUM_ATOMS];
    outp[index1 + 2] = uinp[index1 + 2] * invPoli - fieldScale * fieldp[i + PADDED_NUM_ATOMS * 2];
  }
}

/**
 * This kernel computes [result1 = alpha.source1] and [result2 = alpha.source2].
 *
 * Parameters:
 * float, dimension(*) :: polarity [alpha]
 * real, dimension(3,*) :: source1,source2,result1,result2
 */
extern "C" __global__ void tcgAlpha22(const real* __restrict__ source1, const real* __restrict__ source2,
  const float* __restrict__ polarity, real* __restrict__ result1, real* __restrict__ result2)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < NUM_ATOMS; i += blockDim.x * gridDim.x) {
    int index1 = 3 * i;
    const float poli = polarity[i];
    result1[index1 + 0] = poli * source1[index1 + 0];
    result1[index1 + 1] = poli * source1[index1 + 1];
    result1[index1 + 2] = poli * source1[index1 + 2];
    result2[index1 + 0] = poli * source2[index1 + 0];
    result2[index1 + 1] = poli * source2[index1 + 1];
    result2[index1 + 2] = poli * source2[index1 + 2];
  }
}

/**
 * This is a helper kernel for [output = T.alpha.x] product, where T = 1/alpha + Tu,
 * so [output = x + Tu.(alpha.x)].
 *
 * [-Tu.(alpha.x) = field] is pre-calculated by the induced dipole field kernels,
 * and this kernel will perform [output = x - field].
 *
 * Dimensions:
 * real, dimension(3,*) :: outd,outp [output]
 * real, dimension(3,*) :: xd,xp [x]
 * long long, dimension(PADDED_NUM_ATOMS,3) :: field,fieldp [-Tu.(alpha.x)]
 */
extern "C" __global__ void tcgTAlphaPart3(real* __restrict__ outd, real* __restrict__ outp, const real* __restrict__ xd,
  const real* __restrict__ xp, const long long* __restrict__ field, const long long* __restrict__ fieldp)
{
  const real fieldScale = 1 / (real)0x100000000;
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < NUM_ATOMS; i += blockDim.x * gridDim.x) {
    int index1 = 3 * i;
    outd[index1 + 0] = xd[index1 + 0] - fieldScale * field[i];
    outd[index1 + 1] = xd[index1 + 1] - fieldScale * field[i + PADDED_NUM_ATOMS];
    outd[index1 + 2] = xd[index1 + 2] - fieldScale * field[i + 2 * PADDED_NUM_ATOMS];
    outp[index1 + 0] = xp[index1 + 0] - fieldScale * fieldp[i];
    outp[index1 + 1] = xp[index1 + 1] - fieldScale * fieldp[i + PADDED_NUM_ATOMS];
    outp[index1 + 2] = xp[index1 + 2] - fieldScale * fieldp[i + 2 * PADDED_NUM_ATOMS];
  }
}

/**
 * These macros define wrappers of real* ptr[N] array to pass N pointers to a function by value.
 */
#  define TCG_REAL_POINTERS_FOR_DOT_PRODUCT(NNN) \
    struct TCGRealPtrDT##NNN {                   \
      real* ptr[NNN];                            \
    };
#  define TCG_REAL_POINTERS_FOR_ALPHA_QUADRATIC(NNN) \
    struct TCGRealPtrAQ##NNN {                       \
      real* ptr[NNN];                                \
    };

/**
 * Computes the dot product [scalar = a.b] for N times,
 * and writes the results to outputArray.
 *
 * Parameters:
 * real, dimension(*) :: outputArray
 * RealPointer, dimension(N) :: aPtrArray, bPtrArray
 * Each RealPointer points to a block of dimension(3,*) memory.
 */
#  define TCG_DOT_PRODUCT_N(NNN)                                                                  \
    TCG_REAL_POINTERS_FOR_DOT_PRODUCT(NNN)                                                        \
    extern "C" __global__ void tcgDotProduct##NNN(                                                \
      real* __restrict__ outputArray, TCGRealPtrDT##NNN aPtrArray, TCGRealPtrDT##NNN bPtrArray)   \
    {                                                                                             \
      __shared__ real dotprod_shared__[DOT_PRODUCT_THREAD_BLOCK_SIZE];                            \
      for (int m = 0; m < NNN; ++m) {                                                             \
        if (threadIdx.x == 0) {                                                                   \
          for (int i = 0; i < DOT_PRODUCT_THREAD_BLOCK_SIZE; ++i) {                               \
            dotprod_shared__[i] = 0.0f;                                                           \
          }                                                                                       \
        }                                                                                         \
        __syncthreads();                                                                          \
        const real* __restrict__ a = aPtrArray.ptr[m];                                            \
        const real* __restrict__ b = bPtrArray.ptr[m];                                            \
        for (int index__ = blockIdx.x * blockDim.x + threadIdx.x; index__ < NUM_ATOMS;            \
             index__ += blockDim.x * gridDim.x) {                                                 \
          const int idx__ = 3 * index__;                                                          \
          dotprod_shared__[threadIdx.x]                                                           \
            += (a[idx__] * b[idx__] + a[idx__ + 1] * b[idx__ + 1] + a[idx__ + 2] * b[idx__ + 2]); \
        }                                                                                         \
        __syncthreads();                                                                          \
        if (threadIdx.x == 0) {                                                                   \
          real sum = 0.0f;                                                                        \
          for (int i = 0; i < DOT_PRODUCT_THREAD_BLOCK_SIZE; ++i) {                               \
            sum += dotprod_shared__[i];                                                           \
          }                                                                                       \
          TCG_ATOMIC_ADD_REAL(&outputArray[m], sum);                                              \
        }                                                                                         \
      }                                                                                           \
    }

/**
 * Computes the quadratic form [scalar = a.alpha.b] for N times,
 * and writes the results to outputArray.
 *
 * Parameters:
 * real, dimension(*) :: outputArray
 * float, dimension(*) :: polarity
 * RealPointer, dimension(N) :: aPtrArray, bPtrArray
 * Each RealPointer points to a block of dimension(3,*) memory.
 */
#  define TCG_ALPHA_QUADRATIC_N(NNN)                                                                     \
    TCG_REAL_POINTERS_FOR_ALPHA_QUADRATIC(NNN)                                                           \
    extern "C" __global__ void tcgAlphaQuadratic##NNN(real* __restrict__ outputArray,                    \
      const float* __restrict__ polarity, TCGRealPtrAQ##NNN aPtrArray, TCGRealPtrAQ##NNN bPtrArray)      \
    {                                                                                                    \
      __shared__ real dotprod_shared__[DOT_PRODUCT_THREAD_BLOCK_SIZE];                                   \
      for (int m = 0; m < NNN; ++m) {                                                                    \
        if (threadIdx.x == 0) {                                                                          \
          for (int i = 0; i < DOT_PRODUCT_THREAD_BLOCK_SIZE; ++i) {                                      \
            dotprod_shared__[i] = 0.0f;                                                                  \
          }                                                                                              \
        }                                                                                                \
        __syncthreads();                                                                                 \
        const real* __restrict__ a = aPtrArray.ptr[m];                                                   \
        const real* __restrict__ b = bPtrArray.ptr[m];                                                   \
        for (int index__ = blockIdx.x * blockDim.x + threadIdx.x; index__ < NUM_ATOMS;                   \
             index__ += blockDim.x * gridDim.x) {                                                        \
          const real poli = polarity[index__];                                                           \
          const int idx__ = 3 * index__;                                                                 \
          dotprod_shared__[threadIdx.x]                                                                  \
            += poli * (a[idx__] * b[idx__] + a[idx__ + 1] * b[idx__ + 1] + a[idx__ + 2] * b[idx__ + 2]); \
        }                                                                                                \
        __syncthreads();                                                                                 \
        if (threadIdx.x == 0) {                                                                          \
          real sum = 0.0f;                                                                               \
          for (int i = 0; i < DOT_PRODUCT_THREAD_BLOCK_SIZE; ++i) {                                      \
            sum += dotprod_shared__[i];                                                                  \
          }                                                                                              \
          TCG_ATOMIC_ADD_REAL(&outputArray[m], sum);                                                     \
        }                                                                                                \
      }                                                                                                  \
    }

TCG_DOT_PRODUCT_N(6)
TCG_ALPHA_QUADRATIC_N(6)
TCG_ALPHA_QUADRATIC_N(7)

#  define TCG_SQUARE(x) (x * x)

#  if (!M_TCGGUESS) && M_TCGPREC
/**
 * Parameters:
 * float, dimension(*) :: polarity
 * real, dimension(3,*) :: uind,uinp [for multipole-induced dipole force]
 * real, dimension(3,*) :: uindt,uinpt [for energy]
 * real, dimension(3,PADDED_NUM_ATOMS,tcgnab) :: uad,uap,ubd,ubp [for mutual induced dipole force]
 * real, dimension(3,PADDED_NUM_ATOMS,*) :: workspace,workPolar [pre-computed workspace arrays]
 *
 * real tcgoega = M_TCGOMEGA
 */

extern "C" __global__ void tcgInduce1D1(const float* __restrict__ polarity, real* __restrict__ uind,
  real* __restrict__ uinp, real* __restrict__ uindt, real* __restrict__ uinpt, real* __restrict__ uad,
  real* __restrict__ uap, real* __restrict__ ubd, real* __restrict__ ubp, real* __restrict__ workspace,
  real* __restrict__ workPolar, const real* __restrict__ scalars)
{
#    if M_TCGPEEK && (M_TCGORDER == 1)
  const real tcgomega = M_TCGOMEGA;
#    endif // M_TCGPEEK && (M_TCGORDER == 1)
  // scalars[0 - 5] contain: n0_1, n0_2, t1_1, t1_2, sp0, spp1
  const real n0_1 = scalars[0];
  const real n0_2 = scalars[1];
  const real t1_1 = scalars[2];
  const real t1_2 = scalars[3];
#    if M_TCGORDER == 1
  const real sp0 = scalars[4];
  const real spp1 = scalars[5];
#    endif // M_TCGORDER == 1
  // t4 = n0 / t1
  const real t4_1 = n0_1 / t1_1;
  const real t4_2 = n0_2 / t1_2;
  real* __restrict__ udir = uad;
  real* __restrict__ udirp = uap;
  real* __restrict__ xde_1 = uind;
  real* __restrict__ xde_2 = uinp;
  real* __restrict__ r0_1 = workspace;
  real* __restrict__ r0_2 = workPolar;
  real* __restrict__ P1_1 = workspace + 3 * PADDED_NUM_ATOMS;
  real* __restrict__ P1_2 = workPolar + 3 * PADDED_NUM_ATOMS;
  real* __restrict__ r1_1 = workspace + 3 * PADDED_NUM_ATOMS * 2;
  real* __restrict__ r1_2 = workPolar + 3 * PADDED_NUM_ATOMS * 2;

#    if M_TCGORDER == 1
  real a110_1 = t4_1;
  real a111_1 = 2.0f * sp0 / t1_1;
  real a112_1 = -t4_1 * a111_1;
  real a121_1 = 0.5f * a112_1;

  real a110_2 = t4_2;
  real a111_2 = 2.0f * sp0 / t1_2;
  real a112_2 = -t4_2 * a111_2;
  real a121_2 = 0.5f * a112_2;
#      if M_TCGPEEK
  real a1k10a_1 = tcgomega;
  real a1k11a_1 = -tcgomega * t4_1;
  real a1k11_1 = -2.0f * spp1 * tcgomega / t1_1;
  real a1k12_1 = -t4_1 * a1k11_1;
  real a1k20a_1 = a1k11a_1;
  real a1k21_1 = 0.5f * a1k12_1;

  real a1k10a_2 = tcgomega;
  real a1k11a_2 = -tcgomega * t4_2;
  real a1k11_2 = -2.0f * spp1 * tcgomega / t1_2;
  real a1k12_2 = -t4_2 * a1k11_2;
  real a1k20a_2 = a1k11a_2;
  real a1k21_2 = 0.5f * a1k12_2;
#      endif // M_TCGPEEK
#    endif   // M_TCGORDER == 1

  for (int index = blockIdx.x * blockDim.x + threadIdx.x; index < NUM_ATOMS; index += blockDim.x * gridDim.x) {
    int idx = 3 * index;
    int idx0 = idx;
    int idx1 = idx + 1;
    int idx2 = idx + 2;

    // mu1 = mu0 + gamma0 * p0 (or m0)
    uindt[idx0] += t4_1 * udir[idx0];
    uindt[idx1] += t4_1 * udir[idx1];
    uindt[idx2] += t4_1 * udir[idx2];
    uinpt[idx0] += t4_2 * udirp[idx0];
    uinpt[idx1] += t4_2 * udirp[idx1];
    uinpt[idx2] += t4_2 * udirp[idx2];

#    if (M_TCGORDER > 1) || M_TCGPEEK
    // r1 = r0 - gamma0 * T.p0 (or T.m0)
    r1_1[idx0] = r0_1[idx0] - t4_1 * P1_1[idx0];
    r1_1[idx1] = r0_1[idx1] - t4_1 * P1_1[idx1];
    r1_1[idx2] = r0_1[idx2] - t4_1 * P1_1[idx2];
    r1_2[idx0] = r0_2[idx0] - t4_2 * P1_2[idx0];
    r1_2[idx1] = r0_2[idx1] - t4_2 * P1_2[idx1];
    r1_2[idx2] = r0_2[idx2] - t4_2 * P1_2[idx2];
#    endif // (M_TCGORDER > 1) || M_TCGPEEK

    // TCG1 force and energy
#    if M_TCGORDER == 1
    const float poli = polarity[index];
#      if M_TCGPEEK
    // mu1(peek) = mu1 + omega * alpha * r1
    const float omegaPoli = tcgomega * poli;

    uindt[idx0] += omegaPoli * r1_1[idx0];
    uindt[idx1] += omegaPoli * r1_1[idx1];
    uindt[idx2] += omegaPoli * r1_1[idx2];
    uinpt[idx0] += omegaPoli * r1_2[idx0];
    uinpt[idx1] += omegaPoli * r1_2[idx1];
    uinpt[idx2] += omegaPoli * r1_2[idx2];

    ubp[idx0] = -0.5f * ((a121_1 + a1k21_1) * udir[idx0] + a1k20a_1 * udirp[idx0]);
    ubp[idx1] = -0.5f * ((a121_1 + a1k21_1) * udir[idx1] + a1k20a_1 * udirp[idx1]);
    ubp[idx2] = -0.5f * ((a121_1 + a1k21_1) * udir[idx2] + a1k20a_1 * udirp[idx2]);
    ubd[idx0] = -0.5f * ((a121_2 + a1k21_2) * udirp[idx0] + a1k20a_2 * udir[idx0]);
    ubd[idx1] = -0.5f * ((a121_2 + a1k21_2) * udirp[idx1] + a1k20a_2 * udir[idx1]);
    ubd[idx2] = -0.5f * ((a121_2 + a1k21_2) * udirp[idx2] + a1k20a_2 * udir[idx2]);

    xde_1[idx0] = (a110_2 + a1k10a_2) * r0_1[idx0] + (a111_2 + a1k11_2) * r0_2[idx0] + (a112_2 + a1k12_2) * P1_2[idx0]
      + a1k11a_2 * P1_1[idx0];
    xde_1[idx1] = (a110_2 + a1k10a_2) * r0_1[idx1] + (a111_2 + a1k11_2) * r0_2[idx1] + (a112_2 + a1k12_2) * P1_2[idx1]
      + a1k11a_2 * P1_1[idx1];
    xde_1[idx2] = (a110_2 + a1k10a_2) * r0_1[idx2] + (a111_2 + a1k11_2) * r0_2[idx2] + (a112_2 + a1k12_2) * P1_2[idx2]
      + a1k11a_2 * P1_1[idx2];
    xde_2[idx0] = (a110_1 + a1k10a_1) * r0_2[idx0] + (a111_1 + a1k11_1) * r0_1[idx0] + (a112_1 + a1k12_1) * P1_1[idx0]
      + a1k11a_1 * P1_2[idx0];
    xde_2[idx1] = (a110_1 + a1k10a_1) * r0_2[idx1] + (a111_1 + a1k11_1) * r0_1[idx1] + (a112_1 + a1k12_1) * P1_1[idx1]
      + a1k11a_1 * P1_2[idx1];
    xde_2[idx2] = (a110_1 + a1k10a_1) * r0_2[idx2] + (a111_1 + a1k11_1) * r0_1[idx2] + (a112_1 + a1k12_1) * P1_1[idx2]
      + a1k11a_1 * P1_2[idx2];
#      else  // M_TCGPEEK
    ubp[idx0] = -0.5f * a121_1 * udir[idx0];
    ubp[idx1] = -0.5f * a121_1 * udir[idx1];
    ubp[idx2] = -0.5f * a121_1 * udir[idx2];
    ubd[idx0] = -0.5f * a121_2 * udirp[idx0];
    ubd[idx1] = -0.5f * a121_2 * udirp[idx1];
    ubd[idx2] = -0.5f * a121_2 * udirp[idx2];

    xde_1[idx0] = a110_2 * r0_1[idx0] + a111_2 * r0_2[idx0] + a112_2 * P1_2[idx0];
    xde_1[idx1] = a110_2 * r0_1[idx1] + a111_2 * r0_2[idx1] + a112_2 * P1_2[idx1];
    xde_1[idx2] = a110_2 * r0_1[idx2] + a111_2 * r0_2[idx2] + a112_2 * P1_2[idx2];
    xde_2[idx0] = a110_1 * r0_2[idx0] + a111_1 * r0_1[idx0] + a112_1 * P1_1[idx0];
    xde_2[idx1] = a110_1 * r0_2[idx1] + a111_1 * r0_1[idx1] + a112_1 * P1_1[idx1];
    xde_2[idx2] = a110_1 * r0_2[idx2] + a111_1 * r0_1[idx2] + a112_1 * P1_1[idx2];
#      endif // M_TCGPEEK

    xde_1[idx0] = 0.5f * (poli * xde_1[idx0] + uindt[idx0]);
    xde_1[idx1] = 0.5f * (poli * xde_1[idx1] + uindt[idx1]);
    xde_1[idx2] = 0.5f * (poli * xde_1[idx2] + uindt[idx2]);
    xde_2[idx0] = 0.5f * (poli * xde_2[idx0] + uinpt[idx0]);
    xde_2[idx1] = 0.5f * (poli * xde_2[idx1] + uinpt[idx1]);
    xde_2[idx2] = 0.5f * (poli * xde_2[idx2] + uinpt[idx2]);
#    endif // M_TCGORDER == 1
  }
}

extern "C" __global__ void tcgInduce1D2(const float* __restrict__ polarity, real* __restrict__ uind,
  real* __restrict__ uinp, real* __restrict__ uindt, real* __restrict__ uinpt, real* __restrict__ uad,
  real* __restrict__ uap, real* __restrict__ ubd, real* __restrict__ ubp, real* __restrict__ workspace,
  real* __restrict__ workPolar, const real* __restrict__ scalars)
{
#    if M_TCGPEEK
  const real tcgomega = M_TCGOMEGA;
#    endif // M_TCGPEEK
  // scalars[0 - 5] contain: n0_1, n0_2, t1_1, t1_2, sp0, spp1
  const real n0_1 = scalars[0];
  const real n0_2 = scalars[1];
  const real t1_1 = scalars[2];
  const real t1_2 = scalars[3];
  const real sp0 = scalars[4];
  const real spp1 = scalars[5];
  // t4 = n0 / t1
  const real t4_1 = n0_1 / t1_1;
  const real t4_2 = n0_2 / t1_2;
  real* __restrict__ udir = uad;
  real* __restrict__ udirp = uap;
  real* __restrict__ xde_1 = uind;
  real* __restrict__ xde_2 = uinp;
  real* __restrict__ r0_1 = workspace;
  real* __restrict__ r0_2 = workPolar;
  real* __restrict__ P1_1 = workspace + 3 * PADDED_NUM_ATOMS;
  real* __restrict__ P1_2 = workPolar + 3 * PADDED_NUM_ATOMS;
  real* __restrict__ r1_1 = workspace + 3 * PADDED_NUM_ATOMS * 2;
  real* __restrict__ r1_2 = workPolar + 3 * PADDED_NUM_ATOMS * 2;

  // scalars[6 - 11] contain: n1_1, n1_2, t9_1, t9_2, np1_1, np1_2
  const real n1_1 = scalars[6];
  const real n1_2 = scalars[7];
  const real t9_1 = scalars[8];
  const real t9_2 = scalars[9];
  const real np1_1 = scalars[10];
  const real np1_2 = scalars[11];
  // beta1 = r1.r1/r0.r0 = n1/n0
  // t2 = 1 + beta1
  // t8 = t2*np1 - t4*t9
  // t10 = t1^2 - n0.|P1|^2
  // t3  = t1*t8
  // gamma1 = t10/t3
  const real beta1_1 = n1_1 / n0_1;
  const real beta1_2 = n1_2 / n0_2;
  const real t2_1 = 1.0f + beta1_1;
  const real t2_2 = 1.0f + beta1_2;
  const real t8_1 = t2_1 * np1_1 - t4_1 * t9_1;
  const real t8_2 = t2_2 * np1_2 - t4_2 * t9_2;
  const real t10_1 = TCG_SQUARE(t1_1) - n0_1 * np1_1;
  const real t10_2 = TCG_SQUARE(t1_2) - n0_2 * np1_2;
  const real t3_1 = t1_1 * t8_1;
  const real t3_2 = t1_2 * t8_2;
  const real gamma1_1 = t10_1 / t3_1;
  const real gamma1_2 = t10_2 / t3_2;
  // cross terms
  // sp1 = P1.M.E = P1.udir = spp1
  // b1 = sp0 - gamma1*sp1
  // b2 = sp0*t2 - t4*sp1
  const real sp1 = spp1;
  const real b1_1 = sp0 - gamma1_1 * sp1;
  const real b1_2 = sp0 - gamma1_2 * sp1;
  const real b2_1 = sp0 * t2_1 - t4_1 * sp1;
  const real b2_2 = sp0 * t2_2 - t4_2 * sp1;
#    if M_TCGPEEK
  // scalar[12] contains: spp2
  const real spp2 = scalars[12];
#    endif

  real* __restrict__ uad_2 = uad + 3 * PADDED_NUM_ATOMS;
  real* __restrict__ uap_2 = uap + 3 * PADDED_NUM_ATOMS;
  real* __restrict__ ubd_2 = ubd + 3 * PADDED_NUM_ATOMS;
  real* __restrict__ ubp_2 = ubp + 3 * PADDED_NUM_ATOMS;
  real* __restrict__ t2m0_1 = workspace + 3 * PADDED_NUM_ATOMS * 3;
  real* __restrict__ t2m0_2 = workPolar + 3 * PADDED_NUM_ATOMS * 3;
  real* __restrict__ t3m0_1 = workspace + 3 * PADDED_NUM_ATOMS * 4;
  real* __restrict__ t3m0_2 = workPolar + 3 * PADDED_NUM_ATOMS * 4;

#    if M_TCGORDER == 2
  real a232_1 = t1_1 * t4_1 * gamma1_1 * b2_1 / t3_1;
  real a241_1 = a232_1;
  real a231_1 = -n0_1 * b2_1 / t3_1 - 2.0f * t1_1 * t2_1 * gamma1_1 * b2_1 / t3_1 + t4_1 * gamma1_1 * sp0 / t1_1;
  real a223_1 = a232_1;
  real a222_1 = a231_1;
  real a221_1 = -t4_1 * b1_1 / t1_1 + 2.0f * t1_1 * b2_1 / t3_1 - t4_1 * t9_1 * gamma1_1 * b2_1 / t3_1
    + 2.0f * t2_1 * np1_1 * gamma1_1 * b2_1 / t3_1 - t8_1 * gamma1_1 * b2_1 / t3_1
    - 2.0f * t4_1 * np1_1 * sp0 * gamma1_1 / TCG_SQUARE(t1_1);
  real a220_1 = -gamma1_1 * t4_1;
  real a214_1 = 2.0f * a232_1;
  real a213_1 = 2.0f * a231_1;
  real a212_1 = 2.0f * a221_1;
  real a211_1 = 2.0f
    * (b1_1 / t1_1 - np1_1 * b2_1 / t3_1 - TCG_SQUARE(np1_1) * gamma1_1 * b2_1 / t3_1 / t1_1
        + t9_1 * gamma1_1 * b2_1 / t3_1 + np1_1 * sp0 * gamma1_1 / TCG_SQUARE(t1_1));
  real a21n1_1 = a220_1;
  real a210_1 = t4_1 + gamma1_1 * t2_1;

  real a232_2 = t1_2 * t4_2 * gamma1_2 * b2_2 / t3_2;
  real a241_2 = a232_2;
  real a231_2 = -n0_2 * b2_2 / t3_2 - 2.0f * t1_2 * t2_2 * gamma1_2 * b2_2 / t3_2 + t4_2 * gamma1_2 * sp0 / t1_2;
  real a223_2 = a232_2;
  real a222_2 = a231_2;
  real a221_2 = -t4_2 * b1_2 / t1_2 + 2.0f * t1_2 * b2_2 / t3_2 - t4_2 * t9_2 * gamma1_2 * b2_2 / t3_2
    + 2.0f * t2_2 * np1_2 * gamma1_2 * b2_2 / t3_2 - t8_2 * gamma1_2 * b2_2 / t3_2
    - 2.0f * t4_2 * np1_2 * sp0 * gamma1_2 / TCG_SQUARE(t1_2);
  real a220_2 = -gamma1_2 * t4_2;
  real a214_2 = 2.0f * a232_2;
  real a213_2 = 2.0f * a231_2;
  real a212_2 = 2.0f * a221_2;
  real a211_2 = 2.0f
    * (b1_2 / t1_2 - np1_2 * b2_2 / t3_2 - TCG_SQUARE(np1_2) * gamma1_2 * b2_2 / t3_2 / t1_2
        + t9_2 * gamma1_2 * b2_2 / t3_2 + np1_2 * sp0 * gamma1_2 / TCG_SQUARE(t1_2));
  real a21n1_2 = a220_2;
  real a210_2 = t4_2 + gamma1_2 * t2_2;
#      if M_TCGPEEK
  real a2kwt2_1 = tcgomega * (t2_1 * spp1 - t4_1 * spp2);
  real a2kwg1_1 = tcgomega * (spp1 - gamma1_1 * spp2);
  real a2k41_1 = -a2kwt2_1 * t1_1 * t4_1 * gamma1_1 / t3_1;
  real a2k32_1 = a2k41_1;
  real a2k31_1
    = -tcgomega * t4_1 * gamma1_1 * spp1 / t1_1 + a2kwt2_1 * (n0_1 / t3_1 + 2.0f * t1_1 * t2_1 * gamma1_1 / t3_1);
  real a2k30a_1 = tcgomega * gamma1_1 * t4_1;
  real a2k23_1 = a2k41_1;
  real a2k22_1 = a2k31_1;
  real a2k21_1 = 2.0f * t4_1 * np1_1 / TCG_SQUARE(t1_1) * tcgomega * gamma1_1 * spp1
    + a2kwt2_1 * (-2.0f * t1_1 + (t4_1 * t9_1 - 2.0f * np1_1 * t2_1 + t8_1) * gamma1_1) / t3_1 + t4_1 * a2kwg1_1 / t1_1;
  real a2k21a_1 = a2k30a_1;
  real a2k20a_1 = -tcgomega * (gamma1_1 * t2_1 + t4_1);
  real a2k14_1 = 2.0f * a2k41_1;
  real a2k13_1 = 2.0f * a2k22_1;
  real a2k12_1 = 2.0f * a2k21_1;
  real a2k11_1 = -np1_1 / TCG_SQUARE(t1_1) * tcgomega * gamma1_1 * spp1
    + a2kwt2_1 * (np1_1 + TCG_SQUARE(np1_1) * gamma1_1 / t1_1 - t9_1 * gamma1_1) / t3_1 - a2kwg1_1 / t1_1;
  a2k11_1 = 2.0f * a2k11_1;
  real a2k12a_1 = a2k30a_1;
  real a2k11a_1 = a2k20a_1;
  real a2k10a_1 = tcgomega;

  real a2kwt2_2 = tcgomega * (t2_2 * spp1 - t4_2 * spp2);
  real a2kwg1_2 = tcgomega * (spp1 - gamma1_2 * spp2);
  real a2k41_2 = -a2kwt2_2 * t1_2 * t4_2 * gamma1_2 / t3_2;
  real a2k32_2 = a2k41_2;
  real a2k31_2
    = -tcgomega * t4_2 * gamma1_2 * spp1 / t1_2 + a2kwt2_2 * (n0_2 / t3_2 + 2.0f * t1_2 * t2_2 * gamma1_2 / t3_2);
  real a2k30a_2 = tcgomega * gamma1_2 * t4_2;
  real a2k23_2 = a2k41_2;
  real a2k22_2 = a2k31_2;
  real a2k21_2 = 2.0f * t4_2 * np1_2 / TCG_SQUARE(t1_2) * tcgomega * gamma1_2 * spp1
    + a2kwt2_2 * (-2.0f * t1_2 + (t4_2 * t9_2 - 2.0f * np1_2 * t2_2 + t8_2) * gamma1_2) / t3_2 + t4_2 * a2kwg1_2 / t1_2;
  real a2k21a_2 = a2k30a_2;
  real a2k20a_2 = -tcgomega * (gamma1_2 * t2_2 + t4_2);
  real a2k14_2 = 2.0f * a2k41_2;
  real a2k13_2 = 2.0f * a2k22_2;
  real a2k12_2 = 2.0f * a2k21_2;
  real a2k11_2 = -np1_2 / TCG_SQUARE(t1_2) * tcgomega * gamma1_2 * spp1
    + a2kwt2_2 * (np1_2 + TCG_SQUARE(np1_2) * gamma1_2 / t1_2 - t9_2 * gamma1_2) / t3_2 - a2kwg1_2 / t1_2;
  a2k11_2 = 2.0f * a2k11_2;
  real a2k12a_2 = a2k30a_2;
  real a2k11a_2 = a2k20a_2;
  real a2k10a_2 = tcgomega;
#      endif // M_TCGPEEK
#    endif   // M_TCGORDER == 2

  for (int index = blockIdx.x * blockDim.x + threadIdx.x; index < NUM_ATOMS; index += blockDim.x * gridDim.x) {
    int idx = 3 * index;
    int idx0 = idx;
    int idx1 = idx + 1;
    int idx2 = idx + 2;

    // mu2 = mu1 + gamma1*p1 = mu1 + gamma1*(M.r1 + beta1*p0)
    //     = mu1 + gamma1*(t2*p0 - t4*M.T.p0)
    //     = mu1 + gamma1*(t2*M.r0 - t4*M.T.M.r0)
    //     = mu1 + gamma1*M.(t2*r0 - t4*P1)
    const float poli = polarity[index];
    const real gpoli_1 = gamma1_1 * poli;
    const real gpoli_2 = gamma1_2 * poli;

    uindt[idx0] += (t2_1 * r0_1[idx0] - t4_1 * P1_1[idx0]) * gpoli_1;
    uindt[idx1] += (t2_1 * r0_1[idx1] - t4_1 * P1_1[idx1]) * gpoli_1;
    uindt[idx2] += (t2_1 * r0_1[idx2] - t4_1 * P1_1[idx2]) * gpoli_1;
    uinpt[idx0] += (t2_2 * r0_2[idx0] - t4_2 * P1_2[idx0]) * gpoli_2;
    uinpt[idx1] += (t2_2 * r0_2[idx1] - t4_2 * P1_2[idx1]) * gpoli_2;
    uinpt[idx2] += (t2_2 * r0_2[idx2] - t4_2 * P1_2[idx2]) * gpoli_2;

#    if (M_TCGORDER > 2) || M_TCGPEEK
    // r2 = r1 - gamma1 * T.p1 = r1 - gamma1 * P2
    //    = r1 - gamma1 * (t2*T.M.r0 - t4*T.M.T.M.r0)
    //    = r1 - gamma1 * (t2*P1 - t4*t2m0)
    // reuse r1 as r2
    r1_1[idx0] -= gamma1_1 * (t2_1 * P1_1[idx0] - t4_1 * t2m0_1[idx0]);
    r1_1[idx1] -= gamma1_1 * (t2_1 * P1_1[idx1] - t4_1 * t2m0_1[idx1]);
    r1_1[idx2] -= gamma1_1 * (t2_1 * P1_1[idx2] - t4_1 * t2m0_1[idx2]);
    r1_2[idx0] -= gamma1_2 * (t2_2 * P1_2[idx0] - t4_2 * t2m0_2[idx0]);
    r1_2[idx1] -= gamma1_2 * (t2_2 * P1_2[idx1] - t4_2 * t2m0_2[idx1]);
    r1_2[idx2] -= gamma1_2 * (t2_2 * P1_2[idx2] - t4_2 * t2m0_2[idx2]);
#    endif // (M_TCGORDER >= 2) || M_TCGPEEK

// TCG2 force and energy
#    if M_TCGORDER == 2
#      if M_TCGPEEK
    // mu2(peek) = mu2 + omega * alpha.r2
    const float omegaPoli = tcgomega * poli;

    uindt[idx0] += omegaPoli * r1_1[idx0];
    uindt[idx1] += omegaPoli * r1_1[idx1];
    uindt[idx2] += omegaPoli * r1_1[idx2];
    uinpt[idx0] += omegaPoli * r1_2[idx0];
    uinpt[idx1] += omegaPoli * r1_2[idx1];
    uinpt[idx2] += omegaPoli * r1_2[idx2];

    ubp[idx0] = -0.5f
      * ((a220_1 + a2k20a_1) * r0_2[idx0] + (a221_1 + a2k21_1) * r0_1[idx0]
          + (a222_1 + a231_1 + a2k22_1 + a2k31_1) * P1_1[idx0] + a2k21a_1 * P1_2[idx0]
          + (a223_1 + a241_1 + a2k23_1 + a2k41_1) * t2m0_1[idx0]);
    ubp[idx1] = -0.5f
      * ((a220_1 + a2k20a_1) * r0_2[idx1] + (a221_1 + a2k21_1) * r0_1[idx1]
          + (a222_1 + a231_1 + a2k22_1 + a2k31_1) * P1_1[idx1] + a2k21a_1 * P1_2[idx1]
          + (a223_1 + a241_1 + a2k23_1 + a2k41_1) * t2m0_1[idx1]);
    ubp[idx2] = -0.5f
      * ((a220_1 + a2k20a_1) * r0_2[idx2] + (a221_1 + a2k21_1) * r0_1[idx2]
          + (a222_1 + a231_1 + a2k22_1 + a2k31_1) * P1_1[idx2] + a2k21a_1 * P1_2[idx2]
          + (a223_1 + a241_1 + a2k23_1 + a2k41_1) * t2m0_1[idx2]);
    ubd[idx0] = -0.5f
      * ((a220_2 + a2k20a_2) * r0_1[idx0] + (a221_2 + a2k21_2) * r0_2[idx0]
          + (a222_2 + a231_2 + a2k22_2 + a2k31_2) * P1_2[idx0] + a2k21a_2 * P1_1[idx0]
          + (a223_2 + a241_2 + a2k23_2 + a2k41_2) * t2m0_2[idx0]);
    ubd[idx1] = -0.5f
      * ((a220_2 + a2k20a_2) * r0_1[idx1] + (a221_2 + a2k21_2) * r0_2[idx1]
          + (a222_2 + a231_2 + a2k22_2 + a2k31_2) * P1_2[idx1] + a2k21a_2 * P1_1[idx1]
          + (a223_2 + a241_2 + a2k23_2 + a2k41_2) * t2m0_2[idx1]);
    ubd[idx2] = -0.5f
      * ((a220_2 + a2k20a_2) * r0_1[idx2] + (a221_2 + a2k21_2) * r0_2[idx2]
          + (a222_2 + a231_2 + a2k22_2 + a2k31_2) * P1_2[idx2] + a2k21a_2 * P1_1[idx2]
          + (a223_2 + a241_2 + a2k23_2 + a2k41_2) * t2m0_2[idx2]);

    ubp_2[idx0] = -0.5f * ((a232_1 + a2k32_1) * P1_1[idx0] + a2k30a_1 * r0_2[idx0]);
    ubp_2[idx1] = -0.5f * ((a232_1 + a2k32_1) * P1_1[idx1] + a2k30a_1 * r0_2[idx1]);
    ubp_2[idx2] = -0.5f * ((a232_1 + a2k32_1) * P1_1[idx2] + a2k30a_1 * r0_2[idx2]);
    ubd_2[idx0] = -0.5f * ((a232_2 + a2k32_2) * P1_2[idx0] + a2k30a_2 * r0_1[idx0]);
    ubd_2[idx1] = -0.5f * ((a232_2 + a2k32_2) * P1_2[idx1] + a2k30a_2 * r0_1[idx1]);
    ubd_2[idx2] = -0.5f * ((a232_2 + a2k32_2) * P1_2[idx2] + a2k30a_2 * r0_1[idx2]);

    xde_1[idx0] = (a210_2 + a2k10a_2) * r0_1[idx0] + (a211_2 + a2k11_2) * r0_2[idx0] + (a21n1_2 + a2k11a_2) * P1_1[idx0]
      + (a212_2 + a2k12_2) * P1_2[idx0] + a2k12a_2 * t2m0_1[idx0] + (a213_2 + a2k13_2) * t2m0_2[idx0]
      + (a214_2 + a2k14_2) * t3m0_2[idx0];
    xde_1[idx1] = (a210_2 + a2k10a_2) * r0_1[idx1] + (a211_2 + a2k11_2) * r0_2[idx1] + (a21n1_2 + a2k11a_2) * P1_1[idx1]
      + (a212_2 + a2k12_2) * P1_2[idx1] + a2k12a_2 * t2m0_1[idx1] + (a213_2 + a2k13_2) * t2m0_2[idx1]
      + (a214_2 + a2k14_2) * t3m0_2[idx1];
    xde_1[idx2] = (a210_2 + a2k10a_2) * r0_1[idx2] + (a211_2 + a2k11_2) * r0_2[idx2] + (a21n1_2 + a2k11a_2) * P1_1[idx2]
      + (a212_2 + a2k12_2) * P1_2[idx2] + a2k12a_2 * t2m0_1[idx2] + (a213_2 + a2k13_2) * t2m0_2[idx2]
      + (a214_2 + a2k14_2) * t3m0_2[idx2];
    xde_2[idx0] = (a210_1 + a2k10a_1) * r0_2[idx0] + (a211_1 + a2k11_1) * r0_1[idx0] + (a21n1_1 + a2k11a_1) * P1_2[idx0]
      + (a212_1 + a2k12_1) * P1_1[idx0] + a2k12a_1 * t2m0_2[idx0] + (a213_1 + a2k13_1) * t2m0_1[idx0]
      + (a214_1 + a2k14_1) * t3m0_1[idx0];
    xde_2[idx1] = (a210_1 + a2k10a_1) * r0_2[idx1] + (a211_1 + a2k11_1) * r0_1[idx1] + (a21n1_1 + a2k11a_1) * P1_2[idx1]
      + (a212_1 + a2k12_1) * P1_1[idx1] + a2k12a_1 * t2m0_2[idx1] + (a213_1 + a2k13_1) * t2m0_1[idx1]
      + (a214_1 + a2k14_1) * t3m0_1[idx1];
    xde_2[idx2] = (a210_1 + a2k10a_1) * r0_2[idx2] + (a211_1 + a2k11_1) * r0_1[idx2] + (a21n1_1 + a2k11a_1) * P1_2[idx2]
      + (a212_1 + a2k12_1) * P1_1[idx2] + a2k12a_1 * t2m0_2[idx2] + (a213_1 + a2k13_1) * t2m0_1[idx2]
      + (a214_1 + a2k14_1) * t3m0_1[idx2];
#      else  // M_TCGPEEK
    ubp[idx0] = -0.5f
      * (a220_1 * r0_2[idx0] + a221_1 * r0_1[idx0] + (a222_1 + a231_1) * P1_1[idx0] + (a223_1 + a241_1) * t2m0_1[idx0]);
    ubp[idx1] = -0.5f
      * (a220_1 * r0_2[idx1] + a221_1 * r0_1[idx1] + (a222_1 + a231_1) * P1_1[idx1] + (a223_1 + a241_1) * t2m0_1[idx1]);
    ubp[idx2] = -0.5f
      * (a220_1 * r0_2[idx2] + a221_1 * r0_1[idx2] + (a222_1 + a231_1) * P1_1[idx2] + (a223_1 + a241_1) * t2m0_1[idx2]);
    ubd[idx0] = -0.5f
      * (a220_2 * r0_1[idx0] + a221_2 * r0_2[idx0] + (a222_2 + a231_2) * P1_2[idx0] + (a223_2 + a241_2) * t2m0_2[idx0]);
    ubd[idx1] = -0.5f
      * (a220_2 * r0_1[idx1] + a221_2 * r0_2[idx1] + (a222_2 + a231_2) * P1_2[idx1] + (a223_2 + a241_2) * t2m0_2[idx1]);
    ubd[idx2] = -0.5f
      * (a220_2 * r0_1[idx2] + a221_2 * r0_2[idx2] + (a222_2 + a231_2) * P1_2[idx2] + (a223_2 + a241_2) * t2m0_2[idx2]);

    ubp_2[idx0] = -0.5f * a232_1 * P1_1[idx0];
    ubp_2[idx1] = -0.5f * a232_1 * P1_1[idx1];
    ubp_2[idx2] = -0.5f * a232_1 * P1_1[idx2];
    ubd_2[idx0] = -0.5f * a232_2 * P1_2[idx0];
    ubd_2[idx1] = -0.5f * a232_2 * P1_2[idx1];
    ubd_2[idx2] = -0.5f * a232_2 * P1_2[idx2];

    xde_1[idx0] = a210_2 * r0_1[idx0] + a211_2 * r0_2[idx0] + a21n1_2 * P1_1[idx0] + a212_2 * P1_2[idx0]
      + a213_2 * t2m0_2[idx0] + a214_2 * t3m0_2[idx0];
    xde_1[idx1] = a210_2 * r0_1[idx1] + a211_2 * r0_2[idx1] + a21n1_2 * P1_1[idx1] + a212_2 * P1_2[idx1]
      + a213_2 * t2m0_2[idx1] + a214_2 * t3m0_2[idx1];
    xde_1[idx2] = a210_2 * r0_1[idx2] + a211_2 * r0_2[idx2] + a21n1_2 * P1_1[idx2] + a212_2 * P1_2[idx2]
      + a213_2 * t2m0_2[idx2] + a214_2 * t3m0_2[idx2];
    xde_2[idx0] = a210_1 * r0_2[idx0] + a211_1 * r0_1[idx0] + a21n1_1 * P1_2[idx0] + a212_1 * P1_1[idx0]
      + a213_1 * t2m0_1[idx0] + a214_1 * t3m0_1[idx0];
    xde_2[idx1] = a210_1 * r0_2[idx1] + a211_1 * r0_1[idx1] + a21n1_1 * P1_2[idx1] + a212_1 * P1_1[idx1]
      + a213_1 * t2m0_1[idx1] + a214_1 * t3m0_1[idx1];
    xde_2[idx2] = a210_1 * r0_2[idx2] + a211_1 * r0_1[idx2] + a21n1_1 * P1_2[idx2] + a212_1 * P1_1[idx2]
      + a213_1 * t2m0_1[idx2] + a214_1 * t3m0_1[idx2];
#      endif // M_TCGPEEK
    xde_1[idx0] = 0.5f * (poli * xde_1[idx0] + uindt[idx0]);
    xde_1[idx1] = 0.5f * (poli * xde_1[idx1] + uindt[idx1]);
    xde_1[idx2] = 0.5f * (poli * xde_1[idx2] + uindt[idx2]);
    xde_2[idx0] = 0.5f * (poli * xde_2[idx0] + uinpt[idx0]);
    xde_2[idx1] = 0.5f * (poli * xde_2[idx1] + uinpt[idx1]);
    xde_2[idx2] = 0.5f * (poli * xde_2[idx2] + uinpt[idx2]);

    ubd[idx0] *= poli;
    ubd[idx1] *= poli;
    ubd[idx2] *= poli;
    ubd_2[idx0] *= poli;
    ubd_2[idx1] *= poli;
    ubd_2[idx2] *= poli;
    ubp[idx0] *= poli;
    ubp[idx1] *= poli;
    ubp[idx2] *= poli;
    ubp_2[idx0] *= poli;
    ubp_2[idx1] *= poli;
    ubp_2[idx2] *= poli;

    uad_2[idx0] = P1_1[idx0] * poli;
    uad_2[idx1] = P1_1[idx1] * poli;
    uad_2[idx2] = P1_1[idx2] * poli;
    uap_2[idx0] = P1_2[idx0] * poli;
    uap_2[idx1] = P1_2[idx1] * poli;
    uap_2[idx2] = P1_2[idx2] * poli;
#    endif // M_TCGORDER == 2
  }
}
#  endif // (!M_TCGGUESS) && M_TCGPREC

#endif

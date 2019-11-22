// write out coordinate 
extern "C" __global__ void computeCoordinate(real4* __restrict__ posq) {
    for (int atom = blockIdx.x*blockDim.x+threadIdx.x; atom < NUM_ATOMS; atom += gridDim.x*blockDim.x) {
      if (atom % 3 == 0)
      {
        printf("%i O %15.8f%15.8f%15.8f 1 %i %i \n",atom+1, posq[atom].x, posq[atom].y, posq[atom].z,atom+2, atom+3);
        __syncthreads();
      }
      else
      {
        printf("%i H %15.8f%15.8f%15.8f 2 %i \n",atom+1, posq[atom].x, posq[atom].y, posq[atom].z, atom/3+1);
        __syncthreads();
      }

    }
}

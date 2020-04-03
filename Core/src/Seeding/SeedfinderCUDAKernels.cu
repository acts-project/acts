#include "Acts/Seeding/SeedfinderCUDAKernels.cuh"
#include "Acts/Utilities/Platforms/CUDA/CuUtils.cu"
#include "Acts/Seeding/IExperimentCuts.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SeedfinderConfig.hpp"
#include <cuda.h>
#include <cuda_runtime.h>
#include <iostream>

__global__ void cuSearchDoublet(const unsigned char* isBottom,
				const float* rMvec, const float* zMvec,
				const int* nSpB, const float* rBvec, const float* zBvec, 
				const float* deltaRMin,const float*deltaRMax,const float*cotThetaMax, 
				const float* collisionRegionMin, const float* collisionRegionMax,
				unsigned char* isCompatible				
				);

__global__ void cuTransformCoordinates(const unsigned char* isBottom,
				       const int *  nSpM,
				       const float* spMmat,				       
				       const int* nSpB,
				       const float* spBmat,
				       float* circBmat);

__global__ void cuSearchTriplet(const int*   offset,
				const int*   nSpM,
				const float* spMmat,
				const int*   nSpB, const float* spBmat,
				const int*   nSpT, const float* spTmat,
				const float* circBmat,
				const float* circTmat,
				const float* maxScatteringAngle2, const float* sigmaScattering,
				const float* minHelixDiameter2, const float* pT2perRadius,
				const float* impactMax,
				const int*   nTopPassLimit,
				int* nTopPass, int* tIndex,
				float* curvatures,
				float* impactparameters				
				);

namespace Acts{

  
  void SeedfinderCUDAKernels::searchDoublet(
			        dim3 grid, dim3 block,
				const unsigned char* isBottom,
				const float* rMvec, const float* zMvec,
				const int* nSpB, const float* rBvec, const float* zBvec, 
				const float* deltaRMin,const float*deltaRMax,const float*cotThetaMax, 
				const float* collisionRegionMin, const float* collisionRegionMax,
				unsigned char* isCompatible  ){
    
  cuSearchDoublet<<< grid, block >>>(isBottom,
				     rMvec, zMvec,
				     nSpB, rBvec, zBvec, 
				     deltaRMin, deltaRMax, cotThetaMax, 
				     collisionRegionMin, collisionRegionMax,
				     isCompatible );
  gpuErrChk( cudaGetLastError() );
  }

  void SeedfinderCUDAKernels::transformCoordinates(
				   dim3 grid, dim3 block,
				   const unsigned char* isBottom,
				   const int *  nSpM,
				   const float* spMmat,
				   const int*   nSpB,
				   const float* spBmat,
				   float* circBmat){
    
    cuTransformCoordinates<<< grid, block >>>(isBottom, nSpM, spMmat, nSpB, spBmat, circBmat);
    gpuErrChk( cudaGetLastError() );  
  }

  void SeedfinderCUDAKernels::searchTriplet(
				dim3 grid, dim3 block,
				const int*   offset,
				const int*   nSpM,
				const float* spMmat,
				const int*   nSpB, const float* spBmat,
				const int*   nSpT, const float* spTmat,
				const float* circBmat,
				const float* circTmat,
				const float* maxScatteringAngle2, const float* sigmaScattering,
				const float* minHelixDiameter2,   const float* pT2perRadius,
				const float* impactMax,           const int*   nTopPassLimit,	  
				int*   nTopPass,   int* tIndex,
				float* curvatures, float* impactparameters
				){
cuSearchTriplet<<< grid, block,
      (sizeof(unsigned char)+2*sizeof(float))*block.x >>>(
			       offset,
			       nSpM,
			       spMmat,
			       nSpB, spBmat,
			       nSpT, spTmat,				     
			       circBmat,circTmat,
			       maxScatteringAngle2, sigmaScattering,
			       minHelixDiameter2, pT2perRadius,
			       impactMax, nTopPassLimit,
			       nTopPass, tIndex,
			       curvatures, impactparameters
			       );
  gpuErrChk( cudaGetLastError() );
  }
  
}

__global__ void cuSearchDoublet(const unsigned char* isBottom,
				const float* rMvec, const float* zMvec,
				const int* nSpB, const float* rBvec, const float* zBvec, 	   
				const float* deltaRMin,const float*deltaRMax,const float*cotThetaMax, 
				const float* collisionRegionMin, const float* collisionRegionMax,
				unsigned char* isCompatible 				
				){
  
  int globalId = threadIdx.x+(*nSpB)*blockIdx.x;
  
  float rB = rBvec[threadIdx.x];
  float zB = zBvec[threadIdx.x];
  float rM = rMvec[blockIdx.x];
  float zM = zMvec[blockIdx.x];  
  
  // Doublet search for bottom hits
  isCompatible[globalId] = true;
  
  if (*isBottom == true){    
    float deltaR = rM - rB;
    if (deltaR > *deltaRMax){
      isCompatible[globalId] = false;
    }
        
    if (deltaR < *deltaRMin){
      isCompatible[globalId] = false;
    }
    
    float cotTheta = (zM - zB)/deltaR;
    if (fabs(cotTheta) > *cotThetaMax){
      isCompatible[globalId] = false;
    }

    float zOrigin = zM - rM*cotTheta;
    if (zOrigin < *collisionRegionMin || zOrigin > *collisionRegionMax){
      isCompatible[globalId] = false;
    }    
  }
  
  // Doublet search for top hits
  else if (*isBottom == false){    
    float deltaR = rB - rM;
    if (deltaR < *deltaRMin){
      isCompatible[globalId] = false;
    }

    if (deltaR > *deltaRMax){
      isCompatible[globalId] = false;
    }
    
    if (isCompatible[globalId] == true){
      float cotTheta = (zB - zM)/deltaR;
      if (fabs(cotTheta) > *cotThetaMax){
	isCompatible[globalId] = false;
      }
      
      float zOrigin = zM - rM*cotTheta;
      if (zOrigin < *collisionRegionMin || zOrigin > *collisionRegionMax){
	isCompatible[globalId] = false;
      }
    }    
  }
  
}


__global__ void cuTransformCoordinates(const unsigned char* isBottom,
				       const int *  nSpM,
				       const float* spMmat,
				       const int*   nSpB,
				       const float* spBmat,
				       float*       circBmat){

  int globalId = threadIdx.x+blockDim.x*blockIdx.x;
  if (globalId>=*nSpB) return;
  
  float xM = spMmat[(*nSpM)*0];
  float yM = spMmat[(*nSpM)*1];
  float zM = spMmat[(*nSpM)*2];
  float rM = spMmat[(*nSpM)*3];
  float varianceRM = spMmat[(*nSpM)*4];
  float varianceZM = spMmat[(*nSpM)*5];
  
  float cosPhiM = xM / rM;
  float sinPhiM = yM / rM;
    
  float xB = spBmat[globalId+(*nSpB)*0];
  float yB = spBmat[globalId+(*nSpB)*1];
  float zB = spBmat[globalId+(*nSpB)*2];
  //float rB = spBmat[globalId+(*nSpB)*3];
  float varianceRB = spBmat[globalId+(*nSpB)*4];
  float varianceZB = spBmat[globalId+(*nSpB)*5];

  float deltaX = xB - xM;
  float deltaY = yB - yM;
  float deltaZ = zB - zM;
  
  // calculate projection fraction of spM->sp vector pointing in same
  // direction as
  // vector origin->spM (x) and projection fraction of spM->sp vector pointing
  // orthogonal to origin->spM (y)
  float x = deltaX * cosPhiM + deltaY * sinPhiM;
  float y = deltaY * cosPhiM - deltaX * sinPhiM;
  // 1/(length of M -> SP)
  float iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
  float iDeltaR = std::sqrt(iDeltaR2);

  int bottomFactor = 1 * (int(!(*isBottom))) - 1 * (int(*isBottom));
  // cot_theta = (deltaZ/deltaR)
  float cot_theta = deltaZ * iDeltaR * bottomFactor;
  // VERY frequent (SP^3) access

  // location on z-axis of this SP-duplet
  float Zo = zM - rM * cot_theta;
  
  // transformation of circle equation (x,y) into linear equation (u,v)
  // x^2 + y^2 - 2x_0*x - 2y_0*y = 0
  // is transformed into
  // 1 - 2x_0*u - 2y_0*v = 0
  // using the following m_U and m_V
  // (u = A + B*v); A and B are created later on  
  float U  = x*iDeltaR2;
  float V  = y*iDeltaR2;
  // error term for sp-pair without correlation of middle space point  
  float Er = ((varianceZM + varianceZB) +
	      (cot_theta * cot_theta) * (varianceRM + varianceRB)) * iDeltaR2;  
  
  circBmat[globalId+(*nSpB)*0] = Zo;
  circBmat[globalId+(*nSpB)*1] = cot_theta;
  circBmat[globalId+(*nSpB)*2] = iDeltaR;
  circBmat[globalId+(*nSpB)*3] = Er;
  circBmat[globalId+(*nSpB)*4] = U;
  circBmat[globalId+(*nSpB)*5] = V; 
}

__global__ void cuSearchTriplet(const int*   offset,
				const int*   nSpM,
				const float* spMmat,
				const int*   nSpB, const float* spBmat,
				const int*   nSpT, const float* spTmat,
				const float* circBmat,
				const float* circTmat,
				const float* maxScatteringAngle2, const float* sigmaScattering,
				const float* minHelixDiameter2,    const float* pT2perRadius,
				const float* impactMax,
				const int*   nTopPassLimit,
				int* nTopPass,
				int* tIndex,
				float* curvatures,
				float* impactparameters
				){
  extern __shared__ float shared[];
  
  float* impact   = (float*)shared;
  float* invHelix = (float*)&impact[blockDim.x];
  unsigned char* isPassed = (unsigned char*)&invHelix[blockDim.x];  
  
  float rM         = spMmat[(*nSpM)*3];
  float varianceRM = spMmat[(*nSpM)*4];
  float varianceZM = spMmat[(*nSpM)*5];

  // Zob values from CPU and CUDA are slightly different
  //float Zob        = circBmat[blockId+(*nSpB)*0];
  float cotThetaB  = circBmat[blockIdx.x+(*nSpB)*1];
  float iDeltaRB   = circBmat[blockIdx.x+(*nSpB)*2];
  float ErB        = circBmat[blockIdx.x+(*nSpB)*3];
  float Ub         = circBmat[blockIdx.x+(*nSpB)*4];
  float Vb         = circBmat[blockIdx.x+(*nSpB)*5];
  float iSinTheta2 = (1. + cotThetaB * cotThetaB);
  float scatteringInRegion2 = (*maxScatteringAngle2) * iSinTheta2;
  scatteringInRegion2 *= (*sigmaScattering) * (*sigmaScattering);

  //float Zot        = circTmat[threadId+(*nSpT)*0];
  float cotThetaT  = circTmat[threadIdx.x+(*nSpT)*1];
  float iDeltaRT   = circTmat[threadIdx.x+(*nSpT)*2];
  float ErT        = circTmat[threadIdx.x+(*nSpT)*3];
  float Ut         = circTmat[threadIdx.x+(*nSpT)*4];
  float Vt         = circTmat[threadIdx.x+(*nSpT)*5];

  // add errors of spB-spM and spM-spT pairs and add the correlation term
  // for errors on spM
  float error2 = ErT + ErB +
    2 * (cotThetaB * cotThetaT * varianceRM + varianceZM) * iDeltaRB * iDeltaRT;
  
  float deltaCotTheta = cotThetaB - cotThetaT;
  float deltaCotTheta2 = deltaCotTheta * deltaCotTheta;
  float error;
  float dCotThetaMinusError2;
  
  isPassed[threadIdx.x] = true;
  
  // if the error is larger than the difference in theta, no need to
  // compare with scattering
  if (deltaCotTheta2 - error2 > 0) {
    deltaCotTheta = fabs(deltaCotTheta);
    // if deltaTheta larger than the scattering for the lower pT cut, skip
    error = sqrt(error2);
    dCotThetaMinusError2 =
      deltaCotTheta2 + error2 - 2 * deltaCotTheta * error;
    // avoid taking root of scatteringInRegion
    // if left side of ">" is positive, both sides of unequality can be
    // squared
    // (scattering is always positive)
    
    if (dCotThetaMinusError2 > scatteringInRegion2) {
      isPassed[threadIdx.x] = false;
    }
  }

  // protects against division by 0
  float dU = Ut - Ub;
  if (dU == 0.) {
    isPassed[threadIdx.x] = false;
  }

  // A and B are evaluated as a function of the circumference parameters
  // x_0 and y_0
  float A = (Vt - Vb) / dU;
  float S2 = 1. + A * A;
  float B = Vb - A * Ub;
  float B2 = B * B;
  // sqrt(S2)/B = 2 * helixradius
  // calculated radius must not be smaller than minimum radius
  if (S2 < B2 * (*minHelixDiameter2)) {
    isPassed[threadIdx.x] = false;
  }
  
  // 1/helixradius: (B/sqrt(S2))/2 (we leave everything squared)
  float iHelixDiameter2 = B2 / S2;
  // calculate scattering for p(T) calculated from seed curvature
  float pT2scatter = 4 * iHelixDiameter2 * (*pT2perRadius);
  // TODO: include upper pT limit for scatter calc
  // convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta)
  // from rad to deltaCotTheta
  float p2scatter = pT2scatter * iSinTheta2;
  // if deltaTheta larger than allowed scattering for calculated pT, skip
  if ((deltaCotTheta2 - error2 > 0) &&
      (dCotThetaMinusError2 >
       p2scatter * (*sigmaScattering) * (*sigmaScattering))) {
    isPassed[threadIdx.x] = false;
  }
  // A and B allow calculation of impact params in U/V plane with linear
  // function
  // (in contrast to having to solve a quadratic function in x/y plane)

  impact[threadIdx.x] = fabs((A - B * rM) * rM);
  invHelix[threadIdx.x] = B / sqrt(S2);

  if (impact[threadIdx.x] > (*impactMax)){
    isPassed[threadIdx.x] = false;
  }    
 
  __syncthreads();
  // The following is an iteration to avoid atomic operation
  // But do not use this as it slows down the kernel execution speed
  /*
  if (threadIdx.x==0){
    for (int i_th=0; i_th<blockDim.x; i_th++){
      if (isPassed[i_th] == true){
	int pos = nTopPass[blockIdx.x];
	impactparameters[pos+(*nTopPassLimit)*blockIdx.x] = impact[i_th];
	curvatures      [pos+(*nTopPassLimit)*blockIdx.x] = invHelix[i_th];
	tIndex          [pos+(*nTopPassLimit)*blockIdx.x] = i_th + (*offset);
	nTopPass[blockIdx.x] +=1;
	if( nTopPass[blockIdx.x] == *nTopPassLimit){
	  break;
	}
      }
    }
  }
  */
  
  // The index will be different (and not deterministic) becuase of atomic operation
  // It will be resorted after kernel call
  if (isPassed[threadIdx.x] == true){
    int pos = atomicAdd(&nTopPass[blockIdx.x],1);
    if (pos<*nTopPassLimit){
      impactparameters[pos+(*nTopPassLimit)*blockIdx.x] = impact[threadIdx.x];
      curvatures      [pos+(*nTopPassLimit)*blockIdx.x] = invHelix[threadIdx.x];
      tIndex          [pos+(*nTopPassLimit)*blockIdx.x] = threadIdx.x + (*offset);      
    }
  }
  
  __syncthreads();
  if (threadIdx.x == 0 && nTopPass[blockIdx.x] > *nTopPassLimit){
    nTopPass[blockIdx.x] = *nTopPassLimit;
  }
  
}

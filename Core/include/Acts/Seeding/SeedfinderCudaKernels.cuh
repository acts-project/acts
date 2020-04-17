// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>

namespace Acts{
  
class SeedfinderCudaKernels {

public: 
  
  static void searchDoublet( const dim3 grid, const dim3 block,
			     //const bool* isLast, const int* offset,
			     const int* nSpM, const float* spMmat,
			     const int* nSpB, const float* spBmat,
			     const int* nSpT, const float* spTmat,
			     const float* deltaRMin,
			     const float* deltaRMax,
			     const float* cotThetaMax, 
			     const float* collisionRegionMin, 
			     const float* collisionRegionMax,
			     int*  nSpMcomp,
			     int*  nSpBcompPerSpM_Max,
			     int*  nSpTcompPerSpM_Max,			     
			     int*  nSpBcompPerSpM,
			     int*  nSpTcompPerSpM,
			     int*  McompIndex,			     
			     int*  BcompIndex,
			     int*  tmpBcompIndex,
			     int*  TcompIndex,
			     int*  tmpTcompIndex
			     );
  
  static void transformCoordinate( const dim3 grid, const dim3 block,
				   const int*   nSpM,
				   const float* spMmat,
				   const int*   McompIndex,
				   const int*   nSpB,
				   const float* spBmat,
				   const int*   nSpBcompPerSpM_Max,
				   const int*   BcompIndex,
				   const int*   nSpT,
				   const float* spTmat,
				   const int*   nSpTcompPerSpM_Max,
				   const int*   TcompIndex,
				   float* spMcompMat,			    
				   float* spBcompMatPerSpM,
				   float* circBcompMatPerSpM,
				   float* spTcompMatPerSpM,
				   float* circTcompMatPerSpM);
  
  static void searchTriplet( const dim3 grid, const dim3 block,
			     const int*   offset,
			     const int*   nSpMcomp,
			     const float* spMcompMat,
			     const int*   nSpBcompPerSpM_Max,
			     const float* circBcompMatPerSpM,
			     const int*   nSpTcompPerSpM_Max,
			     const float* circTcompMatPerSpM,
			     const float* maxScatteringAngle2,
			     const float* sigmaScattering,
			     const float* minHelixDiameter2,
			     const float* pT2perRadius,
			     const float* impactMax,
			     const int*   nTrplPerSpMLimit,
			     int* nTrplPerSpM,
			     int* tIndex,
			     int* bIndex,
			     float* curvatures,
			     float* impactparameters,
			     cudaStream_t* stream
			     );

};
}

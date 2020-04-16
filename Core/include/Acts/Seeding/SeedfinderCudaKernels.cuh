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
			     const bool* isBottom, 
			     const int* nSpM, const float* spMmat,
			     const int* nSpB, const float* spBmat,
			     const float* deltaRMin,
			     const float* deltaRMax,
			     const float* cotThetaMax, 
			     const float* collisionRegionMin, 
			     const float* collisionRegionMax,
			     bool* isCompatible,
			     int*  nSpBcomp
			     );

  static void reduceMatrix( const dim3 grid, const dim3 block,
			    const int*   nSpB,
			    const float* spBmat,
			    const int*   nSpBcompPerSpM_Max,
			    const int*   bIndex,
			    float* spBcompMatPerSpM);
  
  static void transformCoordinates( const dim3 grid, const dim3 block,
				    const bool*  isBottom,
				    const float* spMcompMat,
				    const int*   nSpBcompPerSpM_Max,
				    const float* spBcompMatPerSpM_Max,
				    float* circBcompMatPerSpM);
  
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
			     const int*   nTrplPerSpBLimit,
			     int* nTrplPerSpM,
			     int* nTrplPerSpB,
			     int* tIndex,
			     int* bIndex,
			     float* curvatures,
			     float* impactparameters,
			     cudaStream_t* stream
			     );

};
}

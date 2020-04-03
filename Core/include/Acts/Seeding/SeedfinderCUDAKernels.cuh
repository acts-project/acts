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
class CuSeedfinderConfig;
class CuIExperimentCuts;
  
class SeedfinderCUDAKernels {

public: 
  
  static void searchDoublet( dim3 grid, dim3 block,
			     const unsigned char* isBottom, 
			     const float* rMvec, const float* zMvec,
			     const int* nSpB, const float* rBvec, const float* zBvec, 
			     const float* deltaRMin, const float* deltaRMax,
			     const float* cotThetaMax, 
			     const float* collisionRegionMin, 
			     const float* collisionRegionMax,
			     unsigned char* isCompatible 			     
			     );

  static void transformCoordinates( dim3 grid, dim3 block,
				    const unsigned char*  isBottom,
				    const int *  nSpM,
				    const float* spMmat,
				    const int*   nSpB,
				    const float* spBmat,
				    float* circBmat);
  
  static void searchTriplet( dim3 grid, dim3 block,
			     const int*   offset,
			     const int*   nSpM,
			     const float* spMmat,
			     const int*   nSpB, const float* rBmat,
			     const int*   nSpT, const float* rTmat,			     
			     const float* circBmat,
			     const float* circTmat,
			     const float* maxScatteringAngle2, const float* sigmaScattering,
			     const float* minHelixDiameter2,   const float* pT2perRadius,
			     const float* impactMax,           const int*   nTopPassLimit,
			     int*   nTopPass,   int* tIndex,
			     float* curvatures, float* impactparameters
			     );
  
private:
  
};
}

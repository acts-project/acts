// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Seeding/Types.hpp"
#include "Acts/Plugins/Cuda/Utilities/DeviceMatrix.hpp"
#include "Acts/Plugins/Cuda/Utilities/DeviceVector.hpp"
#include "Acts/Plugins/Cuda/Utilities/ResultScalar.hpp"
#include "Acts/Plugins/Cuda/Utilities/StreamWrapper.hpp"

namespace Acts {
namespace Cuda {
namespace details {

void searchDoublet( int grid, int block,
		              int nSpM,
						  const DeviceMatrix<float>& spMmat,
		              int nSpB,
						  const DeviceMatrix<float>& spBmat,
		              int nSpT,
						  const DeviceMatrix<float>& spTmat,
		              float deltaRMin,
		              float deltaRMax,
		              float cotThetaMax,
		              float collisionRegionMin,
		              float collisionRegionMax,
		              ResultScalar<int>& nSpMcomp,
		              ResultScalar<int>& nSpBcompPerSpM_Max,
		              ResultScalar<int>& nSpTcompPerSpM_Max,
		              DeviceVector<int>& nSpBcompPerSpM,
		              DeviceVector<int>& nSpTcompPerSpM,
		              DeviceVector<int>& McompIndex,
		              DeviceMatrix<int>& BcompIndex,
		              DeviceMatrix<int>& tmpBcompIndex,
		              DeviceMatrix<int>& TcompIndex,
		              DeviceMatrix<int>& tmpTcompIndex );

void transformCoordinate( int grid, int block,
			                 int nSpM,
			                 const DeviceMatrix<float>& spMmat,
			                 const DeviceVector<int>& McompIndex,
			                 int nSpB,
			                 const DeviceMatrix<float>& spBmat,
			                 int nSpBcompPerSpM_Max,
			                 const DeviceMatrix<int>& BcompIndex,
			                 int nSpT,
			                 const DeviceMatrix<float>& spTmat,
			                 int nSpTcompPerSpM_Max,
			                 const DeviceMatrix<int>& TcompIndex,
			                 DeviceMatrix<float>& spMcompMat,
			                 DeviceMatrix<float>& spBcompMatPerSpM,
			                 DeviceMatrix<float>& circBcompMatPerSpM,
			                 DeviceMatrix<float>& spTcompMatPerSpM,
			                 DeviceMatrix<float>& circTcompMatPerSpM );

void searchTriplet( int grid, int block,
		              const int*   nSpTcompPerSpM,
		              int          nSpMcomp,
		              const float* spMcompMat,
		              int          nSpBcompPerSpM_Max,
		              const int*   BcompIndex,
		              const float* circBcompMatPerSpM,
		              int          nSpTcompPerSpM_Max,
		              const int*   TcompIndex,
		              const float* spTcompMatPerSpM,
		              const float* circTcompMatPerSpM,
		              float        maxScatteringAngle2,
		              float        sigmaScattering,
		              float        minHelixDiameter2,
		              float        pT2perRadius,
		              float        impactMax,
		              int          nTrplPerSpMLimit,
		              int          nTrplPerSpBLimit,
		              float        deltaInvHelixDiameter,
		              float        impactWeightFactor,
		              float        deltaRMin,
		              float        compatSeedWeight,
		              int          compatSeedLimit,
		              int*         nTrplPerSpM,
		              Triplet*     TripletsPerSpM,
		              const StreamWrapper& streamWrapper );

} // namespace details
} // namespace Cuda
} // namespace Acts

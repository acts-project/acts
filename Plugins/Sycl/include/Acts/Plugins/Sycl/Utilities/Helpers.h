// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once


// store SeedfinderConfig data in float array, index it with enum values
enum eConfigData : int {
  eDeltaRMin = 0,
  eDeltaRMax = 1,
  eCotThetaMax = 2,
  eCollisionRegionMin = 3,
  eCollisionRegionMax = 4,
  eMaxScatteringAngle2 = 5,
  eSigmaScattering = 6,
  eMinHelixDiameter2 = 7,
  ePT2perRadius = 8,
  eDeltaInvHelixDiameter = 9,
  eImpactWeightFactor = 10,
  eFilterDeltaRMin = 11,
  eCompatSeedWeight = 12,
  eImpactMax = 13,
};

// store limits for algorithm in separate array
enum eMaxData : int {
  eCompatSeedLimit = 0,
  eMaxSeedsPerSpM = 1
};

struct offloadSpacePoint {
    float x;
    float y;
    float z;
    float r;
    float varR;
    float varZ;
};

// Offload parameters required to calculate circle with linear equation.
struct offloadLinEqCircle{
    float zo;
    float cotTheta;
    float iDeltaR;
    float er;
    float u;
    float v;
};
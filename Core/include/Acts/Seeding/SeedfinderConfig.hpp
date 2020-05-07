// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
// forward declaration to avoid cyclic dependence
template <typename T>
class SeedFilter;

template <typename SpacePoint>
struct SeedfinderConfig {
  std::shared_ptr<Acts::SeedFilter<SpacePoint>> seedFilter;

  // Seed Cuts
  // lower cutoff for seeds in MeV
  // FIXME: Acts units
  float minPt = 400.;
  // cot of maximum theta angle
  // equivalent to 2.7 eta (pseudorapidity)
  float cotThetaMax = 7.40627;
  // minimum distance in mm in r between two measurements within one seed
  float deltaRMin = 5;
  // maximum distance in mm in r between two measurements within one seed
  float deltaRMax = 270;

  // FIXME: this is not used yet
  //        float upperPtResolutionPerSeed = 20* Acts::GeV;

  // the delta for inverse helix radius up to which compared seeds
  // are considered to have a compatible radius. delta of inverse radius
  // leads to this value being the cutoff. unit is 1/mm. default value
  // of 0.00003 leads to all helices with radius>33m to be considered compatible

  // impact parameter in mm
  // FIXME: Acts units
  float impactMax = 20.;

  // how many sigmas of scattering angle should be considered?
  float sigmaScattering = 5;

  // for how many seeds can one SpacePoint be the middle SpacePoint?
  int maxSeedsPerSpM = 5;

  // Geometry Settings
  // Detector ROI
  // limiting location of collision region in z
  float collisionRegionMin = -150;
  float collisionRegionMax = +150;
  float phiMin = -M_PI;
  float phiMax = M_PI;
  // limiting location of measurements
  float zMin = -2800;
  float zMax = 2800;
  float rMax = 600;
  // WARNING: if rMin is smaller than impactMax, the bin size will be 2*pi,
  // which will make seeding very slow!
  float rMin = 33;

  // Unit in kiloTesla
  // FIXME: Acts units
  float bFieldInZ = 0.00208;
  // location of beam in x,y plane.
  // used as offset for Space Points
  Acts::Vector2D beamPos{0, 0};

  // average radiation lengths of material on the length of a seed. used for
  // scattering.
  // default is 5%
  // TODO: necessary to make amount of material dependent on detector region?
  float radLengthPerSeed = 0.05;
  // alignment uncertainties, used for uncertainties in the
  // non-measurement-plane of the modules
  // which otherwise would be 0
  // will be added to spacepoint measurement uncertainties (and therefore also
  // multiplied by sigmaError)
  // FIXME: call align1 and align2
  float zAlign = 0;
  float rAlign = 0;
  // used for measurement (+alignment) uncertainties.
  // find seeds within 5sigma error ellipse
  float sigmaError = 5;

  // derived values, set on Seedfinder construction
  float highland = 0;
  float maxScatteringAngle2 = 0;
  float pTPerHelixRadius = 0;
  float minHelixDiameter2 = 0;
  float pT2perRadius = 0;

  // only for Cuda plugin
  int maxBlockSize = 1024;
  int nTrplPerSpBLimit = 100;
  int nAvgTrplPerSpBLimit = 2;
};
}  // namespace Acts

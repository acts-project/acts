// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/AtlasStepper.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"

#include <limits>

#include "PropagationDatasets.hpp"
#include "PropagationTests.hpp"

namespace {

namespace ds = ActsTests::PropagationDatasets;
using namespace Acts::UnitLiterals;

using MagneticField = Acts::ConstantBField;
using AtlasStepper = Acts::AtlasStepper<MagneticField>;
using AtlasPropagator = Acts::Propagator<AtlasStepper>;
using EigenStepper = Acts::EigenStepper<MagneticField>;
using EigenPropagator = Acts::Propagator<EigenStepper>;

// absolute parameter tolerances for position, direction, and absolute momentum
constexpr auto epsPos = 1_um;
constexpr auto epsDir = 0.125_mrad;
constexpr auto epsMom = 1_eV;
// relative covariance tolerance
constexpr auto epsCov = 0.0125;
constexpr bool showDebug = false;

constexpr auto bz = 2_T;

const Acts::GeometryContext geoCtx;
const Acts::MagneticFieldContext magCtx;
const MagneticField magField(Acts::Vector3D::UnitZ() * bz);
const AtlasPropagator atlasPropagator{AtlasStepper(magField)};
const EigenPropagator eigenPropagator{EigenStepper(magField)};

}  // namespace

BOOST_AUTO_TEST_SUITE(PropagationCompareAtlasEigenConstant)

BOOST_DATA_TEST_CASE(
    Forward,
    ds::phi* ds::theta* ds::absMomentum* ds::chargeNonZero* ds::pathLength, phi,
    theta, p, q, s) {
  runFreePropagationComparisonTest(
      atlasPropagator, eigenPropagator, geoCtx, magCtx,
      makeParametersCurvilinear(phi, theta, p, q), s, epsPos, epsDir, epsMom,
      epsCov, showDebug);
}

// True forward/backward tracks do not work with z cylinders
BOOST_DATA_TEST_CASE(ToCylinder,
                     ds::phi* ds::thetaNoForwardBackward* ds::absMomentum*
                         ds::chargeNonZero* ds::pathLength,
                     phi, theta, p, q, smax) {
  // transverse radius of the track
  double rt = std::abs(std::sin(theta) * p / bz);
  // make sure the cylinder is reachables
  double cylinderRadius = std::min(rt, 0.5 * smax);

  const auto initialParams = makeParametersCurvilinear(phi, theta, p, q);
  const auto targetSurface = makeTargetCylinder(cylinderRadius);
  auto [atlasParams, atlasPath] =
      transportToSurface(atlasPropagator, geoCtx, magCtx, initialParams,
                         *targetSurface, smax, showDebug);
  auto [eigenParams, eigenPath] =
      transportToSurface(eigenPropagator, geoCtx, magCtx, initialParams,
                         *targetSurface, smax, showDebug);

  checkParametersConsistency(atlasParams, eigenParams, epsPos, epsDir, epsMom);
  CHECK_CLOSE_ABS(atlasPath, eigenPath, epsPos);
}

BOOST_AUTO_TEST_SUITE_END()

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/HelixStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"

#include <utility>

#include "PropagationDatasets.hpp"
#include "PropagationTests.hpp"

namespace {

namespace ds = ActsTests::PropagationDatasets;

using namespace Acts;
using namespace UnitLiterals;

using MagneticField = ConstantBField;
using HelixPropagator = Propagator<HelixStepper>;
using EigenStepper = EigenStepper<>;
using EigenPropagator = Propagator<EigenStepper>;

// The forward-backward and Ridders tests are blind to a wrong sign of the
// bending, so compare against an independent integrator.
constexpr auto epsPos = 1_um;
constexpr auto epsTime = 1_um;
constexpr auto epsDir = 0.125_mrad;
constexpr auto epsMom = 1_eV;
constexpr auto epsCov = 0.025;

// A step to the straight-line distance of a surface passes a surface
// perpendicular to the track by about |(q/p) T x B|^2 h^3 / 3, which is up to
// about 2 mm for these datasets. The next step goes back, so accept
// intersections that far behind the track. Only for planes and discs: they
// have one intersection, while the target cylinder can have its second
// intersection within 1 cm behind the track.
struct TargetReached : public SurfaceReached {
  TargetReached() : SurfaceReached(-1_cm) {}
};

const auto geoCtx = GeometryContext::dangerouslyDefaultConstruct();
const MagneticFieldContext magCtx;

inline std::pair<HelixPropagator, EigenPropagator> makePropagators(double bz) {
  auto field = std::make_shared<MagneticField>(Vector3(0.0, 0.0, bz));
  return {HelixPropagator(HelixStepper(field)),
          EigenPropagator(EigenStepper(field))};
}

}  // namespace

BOOST_AUTO_TEST_SUITE(PropagationCompareHelixEigenConstant)

BOOST_DATA_TEST_CASE(Forward,
                     ds::phi* ds::thetaWithoutBeam* ds::absMomentum*
                         ds::chargeNonZero* ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  auto [helixPropagator, eigenPropagator] = makePropagators(bz);
  runForwardComparisonTest(
      helixPropagator, eigenPropagator, geoCtx, magCtx,
      makeParametersCurvilinearWithCovariance(phi, theta, p, q), s, epsPos,
      epsTime, epsDir, epsMom, epsCov);
}

BOOST_DATA_TEST_CASE(ToCylinderAlongZ,
                     ds::phi* ds::thetaWithoutBeam* ds::absMomentum*
                         ds::chargeNonZero* ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  auto [helixPropagator, eigenPropagator] = makePropagators(bz);
  runToSurfaceComparisonTest(
      helixPropagator, eigenPropagator, geoCtx, magCtx,
      makeParametersCurvilinearWithCovariance(phi, theta, p, q), s,
      ZCylinderSurfaceBuilder(), epsPos, epsTime, epsDir, epsMom, epsCov);
}

BOOST_DATA_TEST_CASE(
    ToDisc,
    ds::phiWithoutAmbiguity* ds::thetaWithoutBeam* ds::absMomentum*
        ds::chargeNonZero* ds::pathLength* ds::magneticField,
    phi, theta, p, q, s, bz) {
  auto [helixPropagator, eigenPropagator] = makePropagators(bz);
  runToSurfaceComparisonTest(
      helixPropagator, eigenPropagator, geoCtx, magCtx,
      makeParametersCurvilinearWithCovariance(phi, theta, p, q), s,
      DiscSurfaceBuilder(), epsPos, epsTime, epsDir, epsMom, epsCov,
      TargetReached());
}

BOOST_DATA_TEST_CASE(ToPlane,
                     ds::phi* ds::thetaWithoutBeam* ds::absMomentum*
                         ds::chargeNonZero* ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  auto [helixPropagator, eigenPropagator] = makePropagators(bz);
  runToSurfaceComparisonTest(
      helixPropagator, eigenPropagator, geoCtx, magCtx,
      makeParametersCurvilinearWithCovariance(phi, theta, p, q), s,
      PlaneSurfaceBuilder(), epsPos, epsTime, epsDir, epsMom, epsCov,
      TargetReached());
}

BOOST_DATA_TEST_CASE(ToStrawAlongZ,
                     ds::phi* ds::thetaWithoutBeam* ds::absMomentum*
                         ds::chargeNonZero* ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  auto [helixPropagator, eigenPropagator] = makePropagators(bz);
  runToSurfaceComparisonTest(
      helixPropagator, eigenPropagator, geoCtx, magCtx,
      makeParametersCurvilinearWithCovariance(phi, theta, p, q), s,
      ZStrawSurfaceBuilder(), epsPos, epsTime, epsDir, epsMom, epsCov);
}

BOOST_AUTO_TEST_SUITE_END()

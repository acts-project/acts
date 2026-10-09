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
#include "Acts/Propagator/GenericRungeKuttaStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/RiddersStepper.hpp"

#include "PropagationDatasets.hpp"
#include "PropagationTests.hpp"

namespace {

namespace ds = ActsTests::PropagationDatasets;

using namespace Acts;
using namespace UnitLiterals;

using MagneticField = ConstantBField;
using Stepper = GenericRungeKuttaStepper;
using RiddersStepper = Experimental::RiddersStepper<Stepper>;
using TestPropagator = Propagator<Stepper>;
using TestRiddersPropagator = Propagator<RiddersStepper>;

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

inline TestPropagator makePropagator(double bz) {
  auto magField = std::make_shared<MagneticField>(Vector3(0.0, 0.0, bz));
  Stepper stepper(std::move(magField));
  return TestPropagator(std::move(stepper));
}

inline TestRiddersPropagator makeRiddersPropagator(double bz) {
  auto magField = std::make_shared<MagneticField>(Vector3(0.0, 0.0, bz));
  RiddersStepper stepper(std::move(magField));
  return TestRiddersPropagator(std::move(stepper));
}

}  // namespace

BOOST_AUTO_TEST_SUITE(PropagationGenericRungeKuttaConstant)

BOOST_DATA_TEST_CASE(ForwardBackward,
                     ds::phi* ds::thetaWithoutBeam* ds::absMomentum*
                         ds::chargeNonZero* ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  runForwardBackwardTest(makePropagator(bz), geoCtx, magCtx,
                         makeParametersCurvilinear(phi, theta, p, q), s, epsPos,
                         epsTime, epsDir, epsMom);
}

BOOST_DATA_TEST_CASE(ToCylinderAlongZ,
                     ds::phi* ds::thetaWithoutBeam* ds::absMomentum*
                         ds::chargeNonZero* ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  runToSurfaceTest(makePropagator(bz), geoCtx, magCtx,
                   makeParametersCurvilinear(phi, theta, p, q), s,
                   ZCylinderSurfaceBuilder(), 1_um, 1_um, 0.125_mrad, epsMom);
}

BOOST_DATA_TEST_CASE(ToDisc,
                     ds::phi* ds::thetaWithoutBeam* ds::absMomentum*
                         ds::chargeNonZero* ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  runToSurfaceTest(makePropagator(bz), geoCtx, magCtx,
                   makeParametersCurvilinear(phi, theta, p, q), s,
                   DiscSurfaceBuilder(), 1_um, 1_um, 0.125_mrad, epsMom,
                   TargetReached());
}

BOOST_DATA_TEST_CASE(ToPlane,
                     ds::phi* ds::thetaWithoutBeam* ds::absMomentum*
                         ds::chargeNonZero* ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  runToSurfaceTest(makePropagator(bz), geoCtx, magCtx,
                   makeParametersCurvilinear(phi, theta, p, q), s,
                   PlaneSurfaceBuilder(), 1_um, 1_um, 0.125_mrad, epsMom,
                   TargetReached());
}

BOOST_DATA_TEST_CASE(ToStrawAlongZ,
                     ds::phi* ds::thetaWithoutBeam* ds::absMomentum*
                         ds::chargeNonZero* ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  runToSurfaceTest(makePropagator(bz), geoCtx, magCtx,
                   makeParametersCurvilinear(phi, theta, p, q), s,
                   ZStrawSurfaceBuilder(), 1_um, 1_um, 0.125_mrad, epsMom);
}

BOOST_DATA_TEST_CASE(CovarianceCurvilinear,
                     ds::phi* ds::thetaWithoutBeam* ds::absMomentum*
                         ds::chargeNonZero* ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  runForwardComparisonTest(
      makePropagator(bz), makeRiddersPropagator(bz), geoCtx, magCtx,
      makeParametersCurvilinearWithCovariance(phi, theta, p, q), s, epsPos,
      epsTime, epsDir, epsMom, epsCov);
}

BOOST_DATA_TEST_CASE(
    CovarianceToCylinderAlongZ,
    ds::phiWithoutAmbiguity* ds::thetaWithoutBeam* ds::absMomentum*
        ds::chargeNonZero* ds::pathLength* ds::magneticField,
    phi, theta, p, q, s, bz) {
  runToSurfaceComparisonTest(
      makePropagator(bz), makeRiddersPropagator(bz), geoCtx, magCtx,
      makeParametersCurvilinearWithCovariance(phi, theta, p, q), s,
      ZCylinderSurfaceBuilder(), 1_um, 1_um, 0.125_mrad, epsMom, epsCov);
}

BOOST_DATA_TEST_CASE(CovarianceToDisc,
                     ds::phi* ds::thetaWithoutBeam* ds::absMomentum*
                         ds::chargeNonZero* ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  runToSurfaceComparisonTest(
      makePropagator(bz), makeRiddersPropagator(bz), geoCtx, magCtx,
      makeParametersCurvilinearWithCovariance(phi, theta, p, q), s,
      DiscSurfaceBuilder(), 1_um, 1_um, 0.125_mrad, epsMom, epsCov,
      TargetReached());
}

BOOST_DATA_TEST_CASE(CovarianceToPlane,
                     ds::phi* ds::thetaWithoutBeam* ds::absMomentum*
                         ds::chargeNonZero* ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  runToSurfaceComparisonTest(
      makePropagator(bz), makeRiddersPropagator(bz), geoCtx, magCtx,
      makeParametersCurvilinearWithCovariance(phi, theta, p, q), s,
      PlaneSurfaceBuilder(), 1_um, 1_um, 0.125_mrad, epsMom, epsCov,
      TargetReached());
}

BOOST_DATA_TEST_CASE(CovarianceToStrawAlongZ,
                     ds::phi* ds::thetaWithoutBeam* ds::absMomentum*
                         ds::chargeNonZero* ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  // the numerical covariance transport to straw surfaces does not seem to be
  // stable. use a higher tolerance for now.
  runToSurfaceComparisonTest(
      makePropagator(bz), makeRiddersPropagator(bz), geoCtx, magCtx,
      makeParametersCurvilinearWithCovariance(phi, theta, p, q), s,
      ZStrawSurfaceBuilder(), 1_um, 1_um, 0.125_mrad, epsMom, 0.125);
}

BOOST_AUTO_TEST_SUITE_END()

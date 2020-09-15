// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "ActsFatras/Digitization/DigitizationError.hpp"
#include "ActsFatras/Digitization/UncorrelatedHitSmearer.hpp"
#include <system_error>

namespace ActsFatras {

BOOST_AUTO_TEST_SUITE(Digitization)

namespace {
struct TestSetter {
  Acts::Result<std::pair<double, double>> operator()() {
    return Acts::Result<std::pair<double, double>>(
        std::make_pair<double, double>(1., 3.));
  }
};

struct SterileSmearer {
  Acts::Result<std::pair<double, double>> operator()() {
    return Acts::Result<std::pair<double, double>>(
        std::make_pair<double, double>(0., 0.));
  }
};

struct InvalidSmearer {
  Acts::Result<std::pair<double, double>> operator()() {
    return Acts::Result<std::pair<double, double>>(
        ActsFatras::DigitizationError::SmearError);
  }
};

const auto pid = Barcode().setVertexPrimary(12).setParticle(23);
const auto gid =
    Acts::GeometryIdentifier().setVolume(1).setLayer(2).setSensitive(3);

auto rec = std::make_shared<Acts::RectangleBounds>(1000, 1000);
auto tSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
    Acts::Transform3D::Identity(), rec);

}  // namespace

BOOST_AUTO_TEST_CASE(HitSmearing_SurfaceMeasurements) {
  Acts::GeometryContext geoCtx = Acts::GeometryContext();

  UncorrelatedHitSmearer hitSmearer;

  // some hit position
  auto p4 = Hit::Vector4(3, 2, 0, 10.);
  // before/after four-momenta are the same
  auto m4 = Hit::Vector4(1, 2, 1, 4);
  auto hit = Hit(gid, pid, p4, m4, m4, 12u);

  TestSetter tSmearFnc;
  SterileSmearer tSterileFnc;
  InvalidSmearer tInvalidFnc;

  auto oneDim = hitSmearer.createSurfaceMeasurement<Acts::eBoundLoc0>(
      geoCtx, hit, *tSurface.get(), {tSmearFnc});

  BOOST_CHECK(oneDim.ok());

  const auto& semearedOne = oneDim.value();
  const auto& smearedParametersOne = semearedOne.parameters();
  BOOST_CHECK(smearedParametersOne.rows() == 1);
  CHECK_CLOSE_ABS(smearedParametersOne[0], 4., Acts::s_epsilon);
  const auto& smearedCovarianceOne = semearedOne.covariance();
  CHECK_CLOSE_ABS(smearedCovarianceOne(Acts::eBoundLoc0, Acts::eBoundLoc0), 9.,
                  Acts::s_epsilon);

  auto twoDim =
      hitSmearer.createSurfaceMeasurement<Acts::eBoundLoc0, Acts::eBoundLoc1>(
          geoCtx, hit, *tSurface.get(), {tSmearFnc, tSmearFnc});

  BOOST_CHECK(twoDim.ok());
  const auto& semearedTwo = twoDim.value();
  const auto& smearedParametersTwo = semearedTwo.parameters();
  BOOST_CHECK(smearedParametersTwo.rows() == 2);
  CHECK_CLOSE_ABS(smearedParametersTwo[0], 4., Acts::s_epsilon);
  CHECK_CLOSE_ABS(smearedParametersTwo[1], 3., Acts::s_epsilon);

  const auto& smearedCovarianceTwo = semearedTwo.covariance();
  CHECK_CLOSE_ABS(smearedCovarianceTwo(Acts::eBoundLoc0, Acts::eBoundLoc0), 9.,
                  Acts::s_epsilon);
  CHECK_CLOSE_ABS(smearedCovarianceTwo(Acts::eBoundLoc1, Acts::eBoundLoc1), 9.,
                  Acts::s_epsilon);

  // Check smearing of time
  auto locYTime =
      hitSmearer.createSurfaceMeasurement<Acts::eBoundLoc1, Acts::eBoundTime>(
          geoCtx, hit, *tSurface.get(), {tSmearFnc, tSmearFnc});
  BOOST_CHECK(locYTime.ok());
  const auto& smearedLocyTime = locYTime.value();
  auto smearedParSetLocyTime = smearedLocyTime.parameterSet();
  BOOST_CHECK(smearedParSetLocyTime.contains<Acts::eBoundLoc1>());
  BOOST_CHECK(smearedParSetLocyTime.contains<Acts::eBoundTime>());
  CHECK_CLOSE_ABS(smearedLocyTime.get<Acts::eBoundLoc1>(), 3., Acts::s_epsilon);
  CHECK_CLOSE_ABS(smearedLocyTime.get<Acts::eBoundTime>(), 11.,
                  Acts::s_epsilon);

  // Use sterile smearer to check if direction is properly translated
  auto phiTheta =
      hitSmearer.createSurfaceMeasurement<Acts::eBoundPhi, Acts::eBoundTheta>(
          geoCtx, hit, *tSurface.get(), {tSterileFnc, tSterileFnc});
  BOOST_CHECK(phiTheta.ok());
  auto phiThetaParSet = phiTheta.value().parameterSet();
  BOOST_CHECK(phiThetaParSet.contains<Acts::eBoundPhi>());
  BOOST_CHECK(phiThetaParSet.contains<Acts::eBoundTheta>());
  CHECK_CLOSE_ABS(phiTheta.value().get<Acts::eBoundPhi>(),
                  Acts::VectorHelpers::phi(hit.unitDirection()),
                  Acts::s_epsilon);
  CHECK_CLOSE_ABS(phiTheta.value().get<Acts::eBoundTheta>(),
                  Acts::VectorHelpers::theta(hit.unitDirection()),
                  Acts::s_epsilon);

  // Finally check an invalid smearing
  auto invalidHitFirst =
      hitSmearer.createSurfaceMeasurement<Acts::eBoundPhi, Acts::eBoundTheta>(
          geoCtx, hit, *tSurface.get(), {tInvalidFnc, tSterileFnc});
  BOOST_CHECK(not invalidHitFirst.ok());

  auto invalidHitMiddle =
      hitSmearer.createSurfaceMeasurement<Acts::eBoundLoc0, Acts::eBoundLoc1,
                                          Acts::eBoundPhi, Acts::eBoundTheta>(
          geoCtx, hit, *tSurface.get(),
          {tSterileFnc, tSterileFnc, tInvalidFnc, tSterileFnc});
  BOOST_CHECK(not invalidHitMiddle.ok());

  auto invalidHitLast =
      hitSmearer.createSurfaceMeasurement<Acts::eBoundPhi, Acts::eBoundTheta>(
          geoCtx, hit, *tSurface.get(), {tSterileFnc, tInvalidFnc});
  BOOST_CHECK(not invalidHitLast.ok());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsFatras
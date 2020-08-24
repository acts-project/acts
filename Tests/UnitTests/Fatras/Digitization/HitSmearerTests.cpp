// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/AbstractVolume.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "ActsFatras/Digitization/HitSmearer.hpp"
#include "ActsFatras/EventData/Hit.hpp"
#include <system_error>

namespace ActsFatras {

BOOST_AUTO_TEST_SUITE(Digitization)

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
    return Acts::Result<std::pair<double, double>>(std::error_code());
  }
};

template <typename object_t>
struct TestHit {
  std::shared_ptr<const object_t> hRefOjbect = nullptr;

  Acts::Vector3D hPosition = Acts::Vector3D(3., 2., 0.);
  Acts::Vector3D hDirection = Acts::Vector3D(1., 2., 1.).normalized();
  double hTime = 10.;

  const Acts::Vector3D position() const { return hPosition; }
  const Acts::Vector3D unitDirection() const { return hDirection; }
  double time() const { return hTime; }

  bool operator==(const TestHit& /*other*/) const { return true; }

  const Acts::Surface& referenceSurface() const { return (*hRefOjbect.get()); }
};

BOOST_AUTO_TEST_CASE(HitSmearing_SurfaceMeasurements) {
  Acts::GeometryContext geoCtx = Acts::GeometryContext();

  auto trf = std::make_shared<Acts::Transform3D>(Acts::Transform3D::Identity());
  auto rec = std::make_shared<Acts::RectangleBounds>(1000, 1000);
  auto tSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(trf, rec);

  HitSmearer hitSmearer;

  std::vector<TestHit<Acts::Surface>> testHits(10, TestHit<Acts::Surface>());
  for (auto& tHit : testHits) {
    tHit.hRefOjbect = tSurface;
  }

  TestSetter tSmearFnc;
  SterileSmearer tSterileFnc;
  InvalidSmearer tInvalidFnc;

  auto oneDim =
      hitSmearer.createSurfaceMeasurement<TestHit<Acts::Surface>, Acts::eLOC_X>(
          geoCtx, testHits[0], *tSurface.get(), {tSmearFnc});

  BOOST_CHECK(oneDim.ok());

  const auto& semearedOne = oneDim.value();
  const auto& smearedParametersOne = semearedOne.parameters();
  BOOST_CHECK(smearedParametersOne.rows() == 1);
  CHECK_CLOSE_ABS(smearedParametersOne[0], 4., Acts::s_epsilon);
  const auto& smearedCovarianceOne = semearedOne.covariance();
  CHECK_CLOSE_ABS(smearedCovarianceOne(Acts::eBoundLoc0, Acts::eBoundLoc0), 9.,
                  Acts::s_epsilon);

  auto twoDim = hitSmearer.createSurfaceMeasurement<TestHit<Acts::Surface>,
                                                    Acts::eLOC_X, Acts::eLOC_Y>(
      geoCtx, testHits[1], *tSurface.get(), {tSmearFnc, tSmearFnc});

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
      hitSmearer.createSurfaceMeasurement<TestHit<Acts::Surface>, Acts::eLOC_Y,
                                          Acts::eBoundTime>(
          geoCtx, testHits[2], *tSurface.get(), {tSmearFnc, tSmearFnc});
  BOOST_CHECK(locYTime.ok());
  const auto& smearedLocyTime = locYTime.value();
  auto smearedParSetLocyTime = smearedLocyTime.parameterSet();
  BOOST_CHECK(smearedParSetLocyTime.contains<Acts::eLOC_Y>());
  BOOST_CHECK(smearedParSetLocyTime.contains<Acts::eBoundTime>());
  CHECK_CLOSE_ABS(smearedLocyTime.get<Acts::eLOC_Y>(), 3., Acts::s_epsilon);
  CHECK_CLOSE_ABS(smearedLocyTime.get<Acts::eBoundTime>(), 11.,
                  Acts::s_epsilon);

  // Use sterile smearer to check if direction is properly translated
  auto phiTheta =
      hitSmearer.createSurfaceMeasurement<TestHit<Acts::Surface>,
                                          Acts::eBoundPhi, Acts::eBoundTheta>(
          geoCtx, testHits[3], *tSurface.get(), {tSterileFnc, tSterileFnc});
  BOOST_CHECK(phiTheta.ok());
  auto phiThetaParSet = phiTheta.value().parameterSet();
  BOOST_CHECK(phiThetaParSet.contains<Acts::eBoundPhi>());
  BOOST_CHECK(phiThetaParSet.contains<Acts::eBoundTheta>());
  CHECK_CLOSE_ABS(phiTheta.value().get<Acts::eBoundPhi>(),
                  Acts::VectorHelpers::phi(testHits[3].hDirection),
                  Acts::s_epsilon);
  CHECK_CLOSE_ABS(phiTheta.value().get<Acts::eBoundTheta>(),
                  Acts::VectorHelpers::theta(testHits[3].hDirection),
                  Acts::s_epsilon);

  // Finally check an invalid smearing
  auto invalidHitFirst =
      hitSmearer.createSurfaceMeasurement<TestHit<Acts::Surface>,
                                          Acts::eBoundPhi, Acts::eBoundTheta>(
          geoCtx, testHits[4], *tSurface.get(), {tInvalidFnc, tSterileFnc});
  BOOST_CHECK(not invalidHitFirst.ok());

  auto invalidHitMiddle =
      hitSmearer.createSurfaceMeasurement<TestHit<Acts::Surface>,
                                          Acts::eBoundLoc0, Acts::eBoundLoc1,
                                          Acts::eBoundPhi, Acts::eBoundTheta>(
          geoCtx, testHits[6], *tSurface.get(),
          {tSterileFnc, tSterileFnc, tInvalidFnc, tSterileFnc});
  BOOST_CHECK(not invalidHitMiddle.ok());

  auto invalidHitLast =
      hitSmearer.createSurfaceMeasurement<TestHit<Acts::Surface>,
                                          Acts::eBoundPhi, Acts::eBoundTheta>(
          geoCtx, testHits[6], *tSurface.get(), {tSterileFnc, tInvalidFnc});
  BOOST_CHECK(not invalidHitLast.ok());
}

/// This is testing the consistency of the Volume smearing function
///
/// As this performs the exact same templated code, only a single test is done
BOOST_AUTO_TEST_CASE(HitSmearing_VolumeMeasurements) {
  auto trf = std::make_shared<Acts::Transform3D>(Acts::Transform3D::Identity());
  auto rec = std::make_shared<Acts::CuboidVolumeBounds>(1000, 1000, 1000);
  auto volume = std::make_shared<const Acts::Volume>(trf, rec);

  HitSmearer hitSmearer;

  std::vector<TestHit<Acts::Volume>> testHits(10, TestHit<Acts::Volume>());
  for (auto& tHit : testHits) {
    tHit.hRefOjbect = volume;
  }

  TestSetter tSmearFnc;

  auto oneDim =
      hitSmearer
          .createVolumeMeasurement<TestHit<Acts::Volume>, Acts::eFreePos0>(
              testHits[0], volume, {tSmearFnc});

  BOOST_CHECK(oneDim.ok());
  const auto& smearedOneDim = oneDim.value();
  CHECK_CLOSE_ABS(smearedOneDim.get<Acts::eFreePos0>(), 4., Acts::s_epsilon);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsFatras
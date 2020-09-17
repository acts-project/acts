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
struct AddSmearer {
  unsigned int stats = 0;

  Acts::Result<std::pair<double, double>> operator()(double value) {
    ++stats;
    return Acts::Result<std::pair<double, double>>(
        std::make_pair<double, double>(value + 1., 3.));
  }
};

struct SterileSmearer {
  Acts::Result<std::pair<double, double>> operator()(double /*ignored*/) {
    return Acts::Result<std::pair<double, double>>(
        std::make_pair<double, double>(0., 0.));
  }
};

struct InvalidSmearer {
  Acts::Result<std::pair<double, double>> operator()(double /*ignored*/) {
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

BOOST_AUTO_TEST_CASE(BoundParameterSmeering) {
  Acts::GeometryContext geoCtx = Acts::GeometryContext();

  UncorrelatedHitSmearer hitSmearer;

  // some hit position
  auto p4 = Hit::Vector4(3, 2, 0, 10.);
  // before/after four-momenta are the same
  auto m4 = Hit::Vector4(1, 2, 1, 4);
  auto hit = Hit(gid, pid, p4, m4, m4, 12u);

  AddSmearer tAddFnc;
  SterileSmearer tSterileFnc;
  InvalidSmearer tInvalidFnc;

  auto oneDim = hitSmearer.smearedParameterSet<Acts::eBoundLoc0>(
      geoCtx, hit, *tSurface.get(), {tAddFnc});

  BOOST_CHECK(oneDim.ok());

  const auto& smearedOne = oneDim.value();
  BOOST_CHECK(smearedOne.contains<Acts::eBoundLoc0>());

  const auto& sParametersOne = smearedOne.getParameters();
  const auto& sCovarianceOne = smearedOne.getCovariance().value();

  CHECK_CLOSE_ABS(sParametersOne[0], 4., Acts::s_epsilon);
  CHECK_CLOSE_ABS(sCovarianceOne(0, 0), 9., Acts::s_epsilon);

  auto twoDim =
      hitSmearer.smearedParameterSet<Acts::eBoundLoc0, Acts::eBoundLoc1>(
          geoCtx, hit, *tSurface.get(), {tAddFnc, tAddFnc});

  BOOST_CHECK(twoDim.ok());
  const auto& smearedTwo = twoDim.value();
  BOOST_CHECK(smearedTwo.contains<Acts::eBoundLoc0>());
  BOOST_CHECK(smearedTwo.contains<Acts::eBoundLoc1>());

  const auto& sParametersTwo = smearedTwo.getParameters();
  const auto& sCovarianceTwo = smearedTwo.getCovariance().value();

  CHECK_CLOSE_ABS(sParametersTwo[0], 4., Acts::s_epsilon);
  CHECK_CLOSE_ABS(sParametersTwo[1], 3., Acts::s_epsilon);
  CHECK_CLOSE_ABS(sCovarianceTwo(0, 0), 9., Acts::s_epsilon);
  CHECK_CLOSE_ABS(sCovarianceTwo(1, 1), 9., Acts::s_epsilon);

  // Check smearing of time
  auto locYTime =
      hitSmearer.smearedParameterSet<Acts::eBoundLoc1, Acts::eBoundTime>(
          geoCtx, hit, *tSurface.get(), {tAddFnc, tAddFnc});
  BOOST_CHECK(locYTime.ok());
  const auto& smearedLocyTime = locYTime.value();
  BOOST_CHECK(smearedLocyTime.contains<Acts::eBoundLoc1>());
  BOOST_CHECK(smearedLocyTime.contains<Acts::eBoundTime>());
  CHECK_CLOSE_ABS(smearedLocyTime.getParameter<Acts::eBoundLoc1>(), 3.,
                  Acts::s_epsilon);
  CHECK_CLOSE_ABS(smearedLocyTime.getParameter<Acts::eBoundTime>(), 11.,
                  Acts::s_epsilon);

  // Use sterile smearer to check if direction is properly translated
  auto phiTheta =
      hitSmearer.smearedParameterSet<Acts::eBoundPhi, Acts::eBoundTheta>(
          geoCtx, hit, *tSurface.get(), {tSterileFnc, tSterileFnc});
  BOOST_CHECK(phiTheta.ok());
  auto phiThetaParSet = phiTheta.value();
  BOOST_CHECK(phiThetaParSet.contains<Acts::eBoundPhi>());
  BOOST_CHECK(phiThetaParSet.contains<Acts::eBoundTheta>());
  CHECK_CLOSE_ABS(phiThetaParSet.getParameter<Acts::eBoundPhi>(),
                  Acts::VectorHelpers::phi(hit.unitDirection()),
                  Acts::s_epsilon);
  CHECK_CLOSE_ABS(phiThetaParSet.getParameter<Acts::eBoundTheta>(),
                  Acts::VectorHelpers::theta(hit.unitDirection()),
                  Acts::s_epsilon);

  // Finally check an invalid smearing
  auto invalidHitFirst =
      hitSmearer.smearedParameterSet<Acts::eBoundPhi, Acts::eBoundTheta>(
          geoCtx, hit, *tSurface.get(), {tInvalidFnc, tSterileFnc});
  BOOST_CHECK(not invalidHitFirst.ok());

  auto invalidHitMiddle =
      hitSmearer.smearedParameterSet<Acts::eBoundLoc0, Acts::eBoundLoc1,
                                     Acts::eBoundPhi, Acts::eBoundTheta>(
          geoCtx, hit, *tSurface.get(),
          {tSterileFnc, tSterileFnc, tInvalidFnc, tSterileFnc});
  BOOST_CHECK(not invalidHitMiddle.ok());

  auto invalidHitLast =
      hitSmearer.smearedParameterSet<Acts::eBoundPhi, Acts::eBoundTheta>(
          geoCtx, hit, *tSurface.get(), {tSterileFnc, tInvalidFnc});
  BOOST_CHECK(not invalidHitLast.ok());
}

BOOST_AUTO_TEST_CASE(FreeParameterSmeering) {
  UncorrelatedHitSmearer hitSmearer;

  // some hit position
  auto p4 = Hit::Vector4(3, 2, 0, 10.);
  // before/after four-momenta are the same
  auto m4 = Hit::Vector4(1, 2, 1, 4);
  auto hit = Hit(gid, pid, p4, m4, m4, 12u);

  AddSmearer tAddFnc;
  SterileSmearer tSterileFnc;

  auto freeSmear =
      hitSmearer.smearedParameterSet<Acts::eFreePos0, Acts::eFreePos1,
                                     Acts::eFreeDir2>(
          hit, {tSterileFnc, tAddFnc, tSterileFnc});

  BOOST_CHECK(freeSmear.ok());
  const auto& freeSet = freeSmear.value();
  BOOST_CHECK(freeSet.contains<Acts::eFreePos0>());
  BOOST_CHECK(freeSet.contains<Acts::eFreePos1>());
  CHECK_CLOSE_ABS(freeSet.getParameter<Acts::eFreePos0>(), 4., Acts::s_epsilon);
  CHECK_CLOSE_ABS(freeSet.getParameter<Acts::eFreePos1>(), 2., Acts::s_epsilon);
  CHECK_CLOSE_ABS(freeSet.getParameter<Acts::eFreeDir2>(),
                  m4.segment<3>(0).normalized()[Acts::eFreeDir2],
                  Acts::s_epsilon);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsFatras
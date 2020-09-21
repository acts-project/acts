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
#include <variant>

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
  Acts::Result<std::pair<double, double>> operator()(double value) {
    return Acts::Result<std::pair<double, double>>(
        std::make_pair<double, double>(value + 0., 0.));
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

  // some hit position
  auto p4 = Hit::Vector4(3, 2, 0, 10.);
  // before/after four-momenta are the same
  auto m4 = Hit::Vector4(1, 2, 1, 4);
  auto hit = Hit(gid, pid, p4, m4, m4, 12u);

  SmearInput sInput{hit, geoCtx, tSurface};

  AddSmearer tAddFnc;
  SterileSmearer tSterileFnc;
  InvalidSmearer tInvalidFnc;

  BoundParametersSmearer<Acts::eBoundLoc0> oneDimSmearer;
  auto oneDim = oneDimSmearer(sInput, {tAddFnc});

  BOOST_CHECK(oneDim.ok());

  const auto& smearedOne = oneDim.value();
  BOOST_CHECK(smearedOne.contains<Acts::eBoundLoc0>());

  const auto& sParametersOne = smearedOne.getParameters();
  const auto& sCovarianceOne = smearedOne.getCovariance().value();

  CHECK_CLOSE_ABS(sParametersOne[0], 4., Acts::s_epsilon);
  CHECK_CLOSE_ABS(sCovarianceOne(0, 0), 9., Acts::s_epsilon);

  BoundParametersSmearer<Acts::eBoundLoc0, Acts::eBoundLoc1> twoDimSmearer;
  auto twoDim = twoDimSmearer(sInput, {tAddFnc, tAddFnc});

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
  BoundParametersSmearer<Acts::eBoundLoc1, Acts::eBoundTime> locYTimeSmearer;
  auto locYTime = locYTimeSmearer(sInput, {tAddFnc, tAddFnc});
  BOOST_CHECK(locYTime.ok());
  const auto& smearedLocyTime = locYTime.value();
  BOOST_CHECK(smearedLocyTime.contains<Acts::eBoundLoc1>());
  BOOST_CHECK(smearedLocyTime.contains<Acts::eBoundTime>());
  CHECK_CLOSE_ABS(smearedLocyTime.getParameter<Acts::eBoundLoc1>(), 3.,
                  Acts::s_epsilon);
  CHECK_CLOSE_ABS(smearedLocyTime.getParameter<Acts::eBoundTime>(), 11.,
                  Acts::s_epsilon);

  // Use sterile BoundParametersSmearer to check if direction is properly
  // translated
  BoundParametersSmearer<Acts::eBoundPhi, Acts::eBoundTheta> phiThetaSmearer;
  auto phiTheta = phiThetaSmearer(sInput, {tSterileFnc, tSterileFnc});
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
  BoundParametersSmearer<Acts::eBoundPhi, Acts::eBoundTheta>
      invalidHitFirstSmearer;
  auto invalidHitFirst =
      invalidHitFirstSmearer(sInput, {tInvalidFnc, tSterileFnc});
  BOOST_CHECK(not invalidHitFirst.ok());

  BoundParametersSmearer<Acts::eBoundLoc0, Acts::eBoundLoc1, Acts::eBoundPhi,
                         Acts::eBoundTheta>
      invalidHitMiddleSmearer;

  auto invalidHitMiddle = invalidHitMiddleSmearer(
      sInput, {tSterileFnc, tSterileFnc, tInvalidFnc, tSterileFnc});
  BOOST_CHECK(not invalidHitMiddle.ok());

  BoundParametersSmearer<Acts::eBoundPhi, Acts::eBoundTheta>
      invalidHitLastSmearer;
  auto invalidHitLast =
      invalidHitLastSmearer(sInput, {tSterileFnc, tInvalidFnc});
  BOOST_CHECK(not invalidHitLast.ok());

  // Test the variant approach for smearers
  using StripSmearer = std::pair<BoundParametersSmearer<Acts::eBoundLoc0>,
                                 std::array<SmearFunction, 1>>;

  using PixelSmearer =
      std::pair<BoundParametersSmearer<Acts::eBoundLoc0, Acts::eBoundLoc1>,
                std::array<SmearFunction, 2>>;

  using SiliconSmearer = std::variant<StripSmearer, PixelSmearer>;

  using SiliconParSets =
      std::variant<Acts::ParameterSet<Acts::BoundIndices, Acts::eBoundLoc0>,
                   Acts::ParameterSet<Acts::BoundIndices, Acts::eBoundLoc0,
                                      Acts::eBoundLoc1>>;

  StripSmearer stripSmearer(oneDimSmearer, {tAddFnc});
  PixelSmearer pixelSmearer(twoDimSmearer, {tAddFnc, tAddFnc});

  std::map<std::string, SiliconSmearer> siSmearers;
  siSmearers["pixel"] = pixelSmearer;
  siSmearers["strip"] = stripSmearer;

  std::map<std::string, Hit> siHits;
  siHits["pixel"] = hit;
  siHits["strip"] = hit;

  std::vector<SiliconParSets> siParSets;

  for (auto& [key, value] : siHits) {
    std::string key_id = key;
    std::visit(
        [&](auto&& sm) {
          auto sParSet = sm.first(sInput, sm.second);
          if (sParSet.ok()) {
            siParSets.push_back(sParSet.value());
            // siParSets.push_Back(sParSet);
            std::cout << "Smearing a hit for " << key_id << std::endl;
            std::cout << " - contains loc0 : "
                      << sParSet.value().template contains<Acts::eBoundLoc0>()
                      << std::endl;
            std::cout << " - contains loc1 : "
                      << sParSet.value().template contains<Acts::eBoundLoc1>()
                      << std::endl;
          }
        },
        siSmearers[key]);
  }

  BOOST_CHECK_EQUAL(siParSets.size(), 2u);
}

BOOST_AUTO_TEST_CASE(FreeParameterSmeering) {
  using FreeSmearer =
      FreeParametersSmearer<Acts::eFreePos0, Acts::eFreePos1, Acts::eFreeDir2>;

  Acts::GeometryContext geoCtx = Acts::GeometryContext();

  FreeSmearer freeSmearer;

  // some hit position
  auto p4 = Hit::Vector4(3, 2, 0, 10.);
  // before/after four-momenta are the same
  auto m4 = Hit::Vector4(1, 2, 1, 4);
  auto hit = Hit(gid, pid, p4, m4, m4, 12u);

  AddSmearer tAddFnc;
  SterileSmearer tSterileFnc;

  auto freeSmearing =
      freeSmearer({hit, geoCtx}, {tSterileFnc, tAddFnc, tSterileFnc});

  BOOST_CHECK(freeSmearing.ok());
  const auto& freeSet = freeSmearing.value();
  BOOST_CHECK(freeSet.contains<Acts::eFreePos0>());
  BOOST_CHECK(freeSet.contains<Acts::eFreePos1>());
  CHECK_CLOSE_ABS(freeSet.getParameter<Acts::eFreePos0>(), 3., Acts::s_epsilon);
  CHECK_CLOSE_ABS(freeSet.getParameter<Acts::eFreePos1>(), 3., Acts::s_epsilon);
  CHECK_CLOSE_ABS(freeSet.getParameter<Acts::eFreeDir2>(),
                  m4.segment<3>(0).normalized()[2], Acts::s_epsilon);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsFatras

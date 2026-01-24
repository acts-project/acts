// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/detail/GenerateParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "ActsExamples/Digitization/Smearers.hpp"
#include "ActsExamples/Io/Json/JsonDigitizationConfig.hpp"
#include "ActsFatras/Digitization/UncorrelatedHitSmearer.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <fstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

using namespace Acts;
using namespace ActsFatras;
using namespace ActsExamples;

namespace {
template <typename generator_t>
struct Fixture {
  generator_t rng;
  // identifiers
  GeometryIdentifier gid;
  Barcode pid;
  // geometry information
  std::shared_ptr<Surface> surface;
  GeometryContext geoCtx = GeometryContext::dangerouslyDefaultConstruct();
  // local and global track parameters
  BoundVector boundParams;
  FreeVector freeParams{};
  // hit information
  Hit hit;

  Fixture(std::uint64_t rngSeed, std::shared_ptr<Surface> surf)
      : rng(rngSeed),
        gid(GeometryIdentifier().withVolume(1).withLayer(2).withSensitive(3)),
        pid(Barcode().withVertexPrimary(12).withParticle(23)),
        surface(std::move(surf)) {
    using namespace UnitLiterals;
    using VectorHelpers::makeVector4;

    surface->assignGeometryId(gid);

    // generate random track parameters
    auto [par, cov] = detail::Test::generateBoundParametersCovariance(rng, {});
    boundParams = par;

    freeParams = transformBoundToFreeParameters(*surface, geoCtx, boundParams);

    // construct hit from free parameters
    Vector4 r4;
    r4.segment<3>(ePos0) = freeParams.segment<3>(eFreePos0);
    r4[eTime] = freeParams[eFreeTime];
    // construct 4-momentum vector assuming m=0
    Vector4 p4;
    p4.segment<3>(eMom0) = freeParams.segment<3>(eFreeDir0).normalized();
    p4[eEnergy] = 1;
    p4 *= std::abs(1_e / freeParams[eFreeQOverP]);
    // same 4-momentum before/after hit
    hit = Hit(gid, pid, r4, p4, p4, 13);
  }
};
}  // namespace

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(GaussianSmearing) {
  nlohmann::json djson = nlohmann::json::parse(R"(
    {
    "acts-geometry-hierarchy-map" : {
    "format-version" : 0,
    "value-identifier" : "digitization-configuration"
    },

  "entries"
      : [
        {
           "volume" : 1,
           "value" : {
            "smearing" : [
              {"index" : 0, "mean" : 0.0, "stddev" : 0.05, "type" : "Gauss", "forcePositiveValues" : true}


            ]
          }
        }
      ]
})");
  double radius = 5.;
  double halfZ = 8.;
  Fixture<RandomEngine> f(
      123567, Surface::makeShared<StrawSurface>(
                  Transform3(Translation3(0., 0., 0.)), radius, halfZ));

  // Get the smearing configuration from the json object
  auto digiConfig =
      DigiConfigConverter("digitization-configuration").fromJson(djson);
  BoundParametersSmearer<RandomEngine, 1u> s;

  for (auto& el : digiConfig) {
    for (auto& smearing : el.smearingDigiConfig.params) {
      // check if the forcePositiveValue parameter is successfully parsed
      BOOST_CHECK(smearing.forcePositiveValues);
      std::fill(std::begin(s.indices), std::end(s.indices),
                static_cast<BoundIndices>(smearing.index));
      std::fill(std::begin(s.smearFunctions), std::end(s.smearFunctions),
                smearing.smearFunction);
      std::fill(std::begin(s.forcePositive), std::end(s.forcePositive),
                smearing.forcePositiveValues);
    }
  }

  auto ret = s(f.rng, f.hit, *f.surface, f.geoCtx);

  BOOST_CHECK(ret.ok());
  auto [par, cov] = ret.value();
  for (std::size_t i = 0; i < s.indices.size(); i++) {
    BOOST_TEST_INFO("Comparing smeared measurement "
                    << i << " originating from bound parameter "
                    << s.indices[i]);
    double ref = f.boundParams[s.indices[i]];
    if (s.forcePositive[i]) {
      ref = std::abs(ref);
    }
    CHECK_CLOSE_REL(par[i], ref, 0.15);
  }
}

BOOST_AUTO_TEST_CASE(DigitizationConfigRoundTrip) {
  std::ofstream out;

  // As all SurfaceBounds have the same streaming API only a one is
  // tested here, all others are tests are identical

  DigiComponentsConfig dcf;

  GeometricConfig gdc;

  BinUtility segmentation;
  segmentation += BinUtility(336, -8.4, 8.4, open, AxisDirection::AxisX);
  segmentation += BinUtility(1280, -36, 36, open, AxisDirection::AxisY);

  gdc.segmentation = segmentation;
  gdc.threshold = 0.01;
  gdc.thickness = 0.15;
  gdc.indices = {eBoundLoc0, eBoundLoc1};
  gdc.chargeSmearer = Digitization::Gauss(1.0);

  DigiComponentsConfig dcRef;
  dcRef.geometricDigiConfig = gdc;

  nlohmann::json dcJsonOut(dcRef);
  out.open("DigiComponentsConfig.json");
  out << dcJsonOut.dump(2);
  out.close();

  auto in = std::ifstream("DigiComponentsConfig.json",
                          std::ifstream::in | std::ifstream::binary);
  BOOST_CHECK(in.good());
  nlohmann::json dcJsonIn;
  in >> dcJsonIn;
  in.close();

  DigiComponentsConfig dcTest(dcJsonIn);
  BOOST_CHECK(dcTest.geometricDigiConfig.indices ==
              dcRef.geometricDigiConfig.indices);
  BOOST_CHECK_EQUAL(dcTest.geometricDigiConfig.segmentation.dimensions(),
                    dcRef.geometricDigiConfig.segmentation.dimensions());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

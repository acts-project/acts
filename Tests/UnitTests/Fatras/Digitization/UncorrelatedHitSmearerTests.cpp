// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/EventData/detail/GenerateParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsExamples/Digitization/DigitizationConfigurator.hpp"
#include "ActsExamples/Digitization/Smearers.hpp"
#include "ActsExamples/Io/Json/JsonDigitizationConfig.hpp"
#include "ActsFatras/Digitization/DigitizationError.hpp"
#include "ActsFatras/Digitization/UncorrelatedHitSmearer.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/Hit.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <memory>
#include <ostream>
#include <random>
#include <utility>

#include <nlohmann/json.hpp>

namespace {

namespace bd = boost::unit_test::data;

using RandomGenerator = std::default_random_engine;

struct SterileSmearer {
  Acts::Result<std::pair<double, double>> operator()(double value,
                                                     RandomGenerator& /*rng*/) {
    return Acts::Result<std::pair<double, double>>(
        std::make_pair<double, double>(value + 0., 0.));
  }
};

struct AddSmearer {
  double offset = 1.0;

  Acts::Result<std::pair<double, double>> operator()(double value,
                                                     RandomGenerator& /*rng*/) {
    return Acts::Result<std::pair<double, double>>(
        std::make_pair<double, double>(value + offset, 3.));
  }
};

struct InvalidSmearer {
  Acts::Result<std::pair<double, double>> operator()(double /*ignored*/,
                                                     RandomGenerator& /*rng*/) {
    return Acts::Result<std::pair<double, double>>(
        ActsFatras::DigitizationError::SmearingError);
  }
};

template <typename generator_t>
struct Fixture {
  generator_t rng;
  // identifiers
  Acts::GeometryIdentifier gid;
  ActsFatras::Barcode pid;
  // geometry information
  std::shared_ptr<Acts::Surface> surface;
  Acts::GeometryContext geoCtx;
  // local and global track parameters
  Acts::BoundVector boundParams;
  Acts::FreeVector freeParams;
  // hit information
  ActsFatras::Hit hit;

  Fixture(uint64_t rngSeed, std::shared_ptr<Acts::Surface> surf)
      : rng(rngSeed),
        gid(Acts::GeometryIdentifier().setVolume(1).setLayer(2).setSensitive(
            3)),
        pid(ActsFatras::Barcode().setVertexPrimary(12).setParticle(23)),
        surface(std::move(surf)) {
    using namespace Acts::UnitLiterals;
    using Acts::VectorHelpers::makeVector4;

    surface->assignGeometryId(gid);

    // generate random track parameters
    auto [par, cov] =
        Acts::detail::Test::generateBoundParametersCovariance(rng);
    boundParams = par;

    freeParams =
        Acts::transformBoundToFreeParameters(*surface, geoCtx, boundParams);

    // construct hit from free parameters
    Acts::Vector4 r4;
    r4.segment<3>(Acts::ePos0) = freeParams.segment<3>(Acts::eFreePos0);
    r4[Acts::eTime] = freeParams[Acts::eFreeTime];
    // construct 4-momentum vector assuming m=0
    Acts::Vector4 p4;
    p4.segment<3>(Acts::eMom0) =
        freeParams.segment<3>(Acts::eFreeDir0).normalized();
    p4[Acts::eEnergy] = 1;
    p4 *= std::abs(1_e / freeParams[Acts::eFreeQOverP]);
    // same 4-momentum before/after hit
    hit = ActsFatras::Hit(gid, pid, r4, p4, p4, 13);
  }
};

// track parameter indices to test smearing with. q/p smearing is not supported
// in either case.
const Acts::BoundIndices boundIndices[] = {
    Acts::eBoundLoc0, Acts::eBoundLoc1,  Acts::eBoundTime,
    Acts::eBoundPhi,  Acts::eBoundTheta,
};
const Acts::FreeIndices freeIndices[] = {
    Acts::eFreePos0, Acts::eFreePos1, Acts::eFreePos2, Acts::eFreeTime,
    Acts::eFreeDir0, Acts::eFreeDir1, Acts::eFreeDir2,
};

constexpr auto tol = 128 * std::numeric_limits<double>::epsilon();

}  // namespace

BOOST_AUTO_TEST_SUITE(FatrasUncorrelatedHitSmearer)

BOOST_DATA_TEST_CASE(Bound1, bd::make(boundIndices), index) {
  Fixture<RandomGenerator> f(
      123, Acts::Surface::makeShared<Acts::PlaneSurface>(
               Acts::Transform3(Acts::Translation3(3, 2, 1))));
  ActsFatras::BoundParametersSmearer<RandomGenerator, 1u> s;
  s.indices = {index};

  // smearing does not do anything
  {
    s.smearFunctions.fill(SterileSmearer{});
    auto ret = s(f.rng, f.hit, *f.surface, f.geoCtx);
    BOOST_CHECK(ret.ok());
    auto [par, cov] = ret.value();
    CHECK_CLOSE_REL(par[0], f.boundParams[index], tol);
  }
  // smearing adds something
  {
    s.smearFunctions.fill(AddSmearer{-42.0});
    auto ret = s(f.rng, f.hit, *f.surface, f.geoCtx);
    BOOST_CHECK(ret.ok());
    auto [par, cov] = ret.value();
    CHECK_CLOSE_REL(par[0], f.boundParams[index] - 42.0, tol);
  }
  // smearing fails
  {
    s.smearFunctions.fill(InvalidSmearer{});
    auto ret = s(f.rng, f.hit, *f.surface, f.geoCtx);
    BOOST_CHECK(!ret.ok());
    BOOST_CHECK(ret.error());
  }
}

BOOST_AUTO_TEST_CASE(BoundAll) {
  Fixture<RandomGenerator> f(
      12356, Acts::Surface::makeShared<Acts::PlaneSurface>(
                 Acts::Transform3(Acts::Translation3(3, 2, 1))));
  // without q/p
  ActsFatras::BoundParametersSmearer<RandomGenerator, std::size(boundIndices)>
      s;
  std::copy(std::begin(boundIndices), std::end(boundIndices),
            s.indices.begin());

  // smearing does not do anything
  {
    s.smearFunctions.fill(SterileSmearer{});
    auto ret = s(f.rng, f.hit, *f.surface, f.geoCtx);
    BOOST_CHECK(ret.ok());
    auto [par, cov] = ret.value();
    for (std::size_t i = 0; i < s.indices.size(); ++i) {
      BOOST_TEST_INFO("Comparing smeared measurement "
                      << i << " originating from bound parameter "
                      << s.indices[i]);
      CHECK_CLOSE_REL(par[i], f.boundParams[s.indices[i]], tol);
    }
  }
  // smearing adds something
  {
    s.smearFunctions.fill(AddSmearer{-23.0});
    auto ret = s(f.rng, f.hit, *f.surface, f.geoCtx);
    BOOST_CHECK(ret.ok());
    auto [par, cov] = ret.value();
    for (std::size_t i = 0; i < s.indices.size(); ++i) {
      BOOST_TEST_INFO("Comparing smeared measurement "
                      << i << " originating from bound parameter "
                      << s.indices[i]);
      CHECK_CLOSE_REL(par[i], f.boundParams[s.indices[i]] - 23.0, tol);
    }
  }
  // one smearer fails
  {
    s.smearFunctions.fill(SterileSmearer{});
    s.smearFunctions[3] = InvalidSmearer{};
    auto ret = s(f.rng, f.hit, *f.surface, f.geoCtx);
    BOOST_CHECK(!ret.ok());
    BOOST_CHECK(ret.error());
  }
}

BOOST_DATA_TEST_CASE(Free1, bd::make(freeIndices), index) {
  Fixture<RandomGenerator> f(
      1234, Acts::Surface::makeShared<Acts::PlaneSurface>(
                Acts::Transform3(Acts::Translation3(3, 2, 1))));
  ActsFatras::FreeParametersSmearer<RandomGenerator, 1u> s;
  s.indices = {index};

  // smearing does not do anything
  {
    s.smearFunctions.fill(SterileSmearer{});
    auto ret = s(f.rng, f.hit);
    BOOST_CHECK(ret.ok());
    auto [par, cov] = ret.value();
    CHECK_CLOSE_REL(par[0], f.freeParams[index], tol);
  }
  // smearing adds something
  {
    s.smearFunctions.fill(AddSmearer{-42.0});
    auto ret = s(f.rng, f.hit);
    BOOST_CHECK(ret.ok());
    auto [par, cov] = ret.value();
    CHECK_CLOSE_REL(par[0], f.freeParams[index] - 42.0, tol);
  }
  // smearing fails
  {
    s.smearFunctions.fill(InvalidSmearer{});
    auto ret = s(f.rng, f.hit);
    BOOST_CHECK(!ret.ok());
    BOOST_CHECK(ret.error());
  }
}

BOOST_AUTO_TEST_CASE(FreeAll) {
  Fixture<RandomGenerator> f(
      123567, Acts::Surface::makeShared<Acts::PlaneSurface>(
                  Acts::Transform3(Acts::Translation3(3, 2, 1))));
  // without q/p
  ActsFatras::FreeParametersSmearer<RandomGenerator, std::size(freeIndices)> s;
  std::copy(std::begin(freeIndices), std::end(freeIndices), s.indices.begin());

  // smearing does not do anything
  {
    s.smearFunctions.fill(SterileSmearer{});
    auto ret = s(f.rng, f.hit);
    BOOST_CHECK(ret.ok());
    auto [par, cov] = ret.value();
    for (std::size_t i = 0; i < s.indices.size(); ++i) {
      BOOST_TEST_INFO("Comparing smeared measurement "
                      << i << " originating from free parameter "
                      << s.indices[i]);
      CHECK_CLOSE_REL(par[i], f.freeParams[s.indices[i]], tol);
    }
  }
  // smearing adds something
  {
    s.smearFunctions.fill(AddSmearer{42.0});
    auto ret = s(f.rng, f.hit);
    BOOST_CHECK(ret.ok());
    auto [par, cov] = ret.value();
    for (std::size_t i = 0; i < s.indices.size(); ++i) {
      BOOST_TEST_INFO("Comparing smeared measurement "
                      << i << " originating from free parameter "
                      << s.indices[i]);
      CHECK_CLOSE_REL(par[i], f.freeParams[s.indices[i]] + 42.0, tol);
    }
  }
  // one smearer fails
  {
    s.smearFunctions.fill(SterileSmearer{});
    s.smearFunctions[3] = InvalidSmearer{};
    auto ret = s(f.rng, f.hit);
    BOOST_CHECK(!ret.ok());
    BOOST_CHECK(ret.error());
  }
}

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
  Fixture<ActsExamples::RandomEngine> f(
      123567,
      Acts::Surface::makeShared<Acts::StrawSurface>(
          Acts::Transform3(Acts::Translation3(0., 0., 0.)), radius, halfZ));

  // Get the smearing configuration from the json object
  auto digiConfig =
      ActsExamples::DigiConfigConverter("digitization-configuration")
          .fromJson(djson);
  ActsFatras::BoundParametersSmearer<ActsExamples::RandomEngine, 1u> s;

  for (auto& el : digiConfig) {
    for (auto& smearing : el.smearingDigiConfig) {
      // check if the forcePositiveValue parameter is successfully parsed
      BOOST_CHECK(smearing.forcePositiveValues);
      std::fill(std::begin(s.indices), std::end(s.indices),
                static_cast<Acts::BoundIndices>(smearing.index));
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

BOOST_AUTO_TEST_SUITE_END()

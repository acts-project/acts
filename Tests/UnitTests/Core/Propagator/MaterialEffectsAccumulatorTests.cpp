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
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Propagator/detail/MaterialEffectsAccumulator.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"
#include "ActsTests/CommonHelpers/PredefinedMaterials.hpp"

#include <utility>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace {

const ParticleHypothesis particle = ParticleHypothesis::pion();
const double momentum = 1_GeV;
const double qOverP = particle.qOverP(momentum, 1);
const double maxXOverX0Step = 0.01;
const Vector3 direction = Vector3::UnitZ();

detail::MaterialEffectsAccumulator makeAccumulator() {
  detail::MaterialEffectsAccumulator accumulator;
  accumulator.initialize(maxXOverX0Step, particle, momentum);
  return accumulator;
}

}  // namespace

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(PropagatorSuite)

BOOST_AUTO_TEST_CASE(MaterialEffectsAccumulatorStepBack) {
  const Material material = makeSilicon();

  auto overstep = makeAccumulator();
  overstep.accumulate(material, 50_mm, qOverP, qOverP);
  BOOST_REQUIRE_NO_THROW(overstep.accumulate(material, -5_mm, qOverP, qOverP));

  auto direct = makeAccumulator();
  direct.accumulate(material, 45_mm, qOverP, qOverP);

  const auto overstepCov = overstep.computeAdditionalFreeCovariance(direction);
  const auto directCov = direct.computeAdditionalFreeCovariance(direction);
  BOOST_REQUIRE(overstepCov.has_value());
  BOOST_REQUIRE(directCov.has_value());

  for (const auto i : {eFreePos0, eFreeDir0, eFreeQOverP}) {
    BOOST_CHECK_GT((*overstepCov)(i, i), 0.);
  }
  CHECK_CLOSE_REL((*overstepCov)(eFreeDir0, eFreeDir0),
                  (*directCov)(eFreeDir0, eFreeDir0), 1e-3);
  CHECK_CLOSE_REL((*overstepCov)(eFreePos0, eFreePos0),
                  (*directCov)(eFreePos0, eFreePos0), 1e-3);
  CHECK_CLOSE_REL((*overstepCov)(eFreePos0, eFreeDir0),
                  (*directCov)(eFreePos0, eFreeDir0), 1e-3);
  CHECK_CLOSE_REL((*overstepCov)(eFreeQOverP, eFreeQOverP),
                  (*directCov)(eFreeQOverP, eFreeQOverP), 1e-3);
}

BOOST_AUTO_TEST_CASE(MaterialEffectsAccumulatorStepBackBeyondStart) {
  auto accumulator = makeAccumulator();
  accumulator.accumulate(makeSilicon(), 1_mm, qOverP, qOverP);
  BOOST_REQUIRE_NO_THROW(
      accumulator.accumulate(makeSilicon(), -2_mm, qOverP, qOverP));
  BOOST_CHECK(accumulator.isVacuum());
  BOOST_CHECK(!accumulator.computeAdditionalFreeCovariance(direction));
}

void checkSameCovariance(detail::MaterialEffectsAccumulator& overstep,
                         detail::MaterialEffectsAccumulator& direct) {
  const auto overstepCov = overstep.computeAdditionalFreeCovariance(direction);
  const auto directCov = direct.computeAdditionalFreeCovariance(direction);
  BOOST_REQUIRE(overstepCov.has_value());
  BOOST_REQUIRE(directCov.has_value());
  for (const auto [i, j] :
       {std::pair{eFreeDir0, eFreeDir0}, std::pair{eFreePos0, eFreePos0},
        std::pair{eFreePos0, eFreeDir0}, std::pair{eFreeQOverP, eFreeQOverP}}) {
    CHECK_CLOSE_REL((*overstepCov)(i, j), (*directCov)(i, j), 1e-3);
  }
}

BOOST_AUTO_TEST_CASE(MaterialEffectsAccumulatorStepBackDenseOverstep) {
  // Overstep from a light into a dense material, then step back. The step back
  // has to remove the dense material, not the average one.
  auto overstep = makeAccumulator();
  overstep.accumulate(makeLiquidArgon(), 100_mm, qOverP, qOverP);
  overstep.accumulate(makeIron(), 1_mm, qOverP, qOverP);
  overstep.accumulate(makeIron(), -1_mm, qOverP, qOverP);

  auto direct = makeAccumulator();
  direct.accumulate(makeLiquidArgon(), 100_mm, qOverP, qOverP);

  checkSameCovariance(overstep, direct);
}

BOOST_AUTO_TEST_CASE(MaterialEffectsAccumulatorStepBackVacuum) {
  // Stepping back through vacuum must not remove any material.
  auto overstep = makeAccumulator();
  overstep.accumulate(makeSilicon(), 10_mm, qOverP, qOverP);
  overstep.accumulate(Material::Vacuum(), 5_mm, qOverP, qOverP);
  overstep.accumulate(Material::Vacuum(), -5_mm, qOverP, qOverP);

  auto direct = makeAccumulator();
  direct.accumulate(makeSilicon(), 10_mm, qOverP, qOverP);

  checkSameCovariance(overstep, direct);
}

BOOST_AUTO_TEST_CASE(MaterialEffectsAccumulatorStepBackMismatchedMaterial) {
  // The step back claims more matter than was accumulated, e.g. for binned
  // volume material sampled at a different position. The result must stay
  // finite.
  auto accumulator = makeAccumulator();
  accumulator.accumulate(makeLiquidArgon(), 10_mm, qOverP, qOverP);
  BOOST_REQUIRE_NO_THROW(
      accumulator.accumulate(makeIron(), -5_mm, qOverP, qOverP));
  const auto cov = accumulator.computeAdditionalFreeCovariance(direction);
  BOOST_REQUIRE(cov.has_value());
  BOOST_CHECK(cov->allFinite());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

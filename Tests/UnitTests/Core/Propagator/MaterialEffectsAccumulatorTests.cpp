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

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(PropagatorMaterialEffectsAccumulator)

// Used to throw "thickness < 0"
BOOST_AUTO_TEST_CASE(BackwardPropagationRemovesMaterialEffects) {
  const Material material = makeSilicon();
  const ParticleHypothesis particleHypothesis = ParticleHypothesis::pion();
  const double momentum = 1_GeV;
  const double qOverP = particleHypothesis.qOverP(momentum, 1);
  const double maxXOverX0Step = 0.1;
  const double pathLength = 50_mm;
  const Vector3 direction = Vector3::UnitZ();

  detail::MaterialEffectsAccumulator forward;
  forward.initialize(maxXOverX0Step, particleHypothesis, momentum);
  forward.accumulate(material, pathLength, qOverP, qOverP);
  const auto forwardCov = forward.computeAdditionalFreeCovariance(direction);
  BOOST_REQUIRE(forwardCov.has_value());

  detail::MaterialEffectsAccumulator backward;
  backward.initialize(maxXOverX0Step, particleHypothesis, momentum);
  BOOST_REQUIRE_NO_THROW(
      backward.accumulate(material, -pathLength, qOverP, qOverP));
  const auto backwardCov = backward.computeAdditionalFreeCovariance(direction);
  BOOST_REQUIRE(backwardCov.has_value());

  // forward adds, backward removes
  BOOST_CHECK_GT((*forwardCov)(eFreeDir1, eFreeDir1), 0.);
  BOOST_CHECK_GT((*forwardCov)(eFreePos1, eFreePos1), 0.);
  BOOST_CHECK_GT((*forwardCov)(eFreeQOverP, eFreeQOverP), 0.);
  CHECK_CLOSE_ABS(*backwardCov, -*forwardCov, 1e-12);
}

// Vacuum contributes nothing either way
BOOST_AUTO_TEST_CASE(VacuumHasNoEffectInEitherDirection) {
  const ParticleHypothesis particleHypothesis = ParticleHypothesis::pion();
  const double momentum = 1_GeV;
  const double qOverP = particleHypothesis.qOverP(momentum, 1);

  for (const double pathLength : {50_mm, -50_mm}) {
    detail::MaterialEffectsAccumulator accumulator;
    accumulator.initialize(0.1, particleHypothesis, momentum);
    BOOST_REQUIRE_NO_THROW(accumulator.accumulate(Material::Vacuum(), pathLength,
                                                  qOverP, qOverP));
    BOOST_CHECK(accumulator.isVacuum());
    BOOST_CHECK(
        !accumulator.computeAdditionalFreeCovariance(Vector3::UnitZ())
             .has_value());
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

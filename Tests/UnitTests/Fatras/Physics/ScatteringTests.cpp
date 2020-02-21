// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <limits>
#include <random>

#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/Physics/Scattering/GaussianMixture.hpp"
#include "ActsFatras/Physics/Scattering/GeneralMixture.hpp"
#include "ActsFatras/Physics/Scattering/Highland.hpp"
#include "Dataset.hpp"

namespace {
constexpr auto eps = std::numeric_limits<double>::epsilon();

// Common test method that will be instantiated for each scattering model.
template <typename Scattering>
void test(const Scattering& scattering, uint32_t seed,
          const ActsFatras::Particle& before) {
  std::ranlux48 gen(seed);
  ActsFatras::Particle after = before;

  const auto outgoing = scattering(gen, Acts::Test::makePercentSlab(), after);
  // scattering leaves absolute energy/momentum unchanged
  CHECK_CLOSE_REL(after.absMomentum(), before.absMomentum(), eps);
  CHECK_CLOSE_REL(after.energy(), before.energy(), eps);
  // scattering has changed the direction
  BOOST_TEST(before.unitDirection().dot(after.unitDirection()) < 1);
  // scattering creates no new particles
  BOOST_TEST(outgoing.empty());
}
}  // namespace

BOOST_AUTO_TEST_SUITE(FatrasScattering)

BOOST_DATA_TEST_CASE(GeneralMixture, Dataset::parameters, pdg, phi, lambda, p,
                     seed) {
  test(ActsFatras::GeneralMixtureScattering(), seed,
       Dataset::makeParticle(pdg, phi, lambda, p));
}

BOOST_DATA_TEST_CASE(GaussianMixture, Dataset::parameters, pdg, phi, lambda, p,
                     seed) {
  test(ActsFatras::GaussianMixtureScattering(), seed,
       Dataset::makeParticle(pdg, phi, lambda, p));
}

BOOST_DATA_TEST_CASE(Highland, Dataset::parameters, pdg, phi, lambda, p, seed) {
  test(ActsFatras::HighlandScattering(), seed,
       Dataset::makeParticle(pdg, phi, lambda, p));
}

BOOST_AUTO_TEST_SUITE_END()

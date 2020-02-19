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
#include "ActsFatras/Physics/Scattering/Scattering.hpp"
#include "Dataset.hpp"

namespace {
constexpr auto eps = std::numeric_limits<double>::epsilon();

using GeneralMixtureScattering =
    ActsFatras::Scattering<ActsFatras::GeneralMixture>;
using GaussianMixtureScattering =
    ActsFatras::Scattering<ActsFatras::GaussianMixture>;
using HighlandScattering = ActsFatras::Scattering<ActsFatras::Highland>;

// Common test method that will be instantiated for each scattering model.
template <typename Scattering>
void test(const Scattering& scattering, uint32_t seed,
          const ActsFatras::Particle& before) {
  std::ranlux48 gen(seed);
  ActsFatras::Particle after = before;

  const auto outgoing = scattering(gen, Acts::Test::makePercentSlab(), after);
  // scattering leaves absolute energy/momentum unchanged
  CHECK_CLOSE_REL(after.momentum(), before.momentum(), eps);
  CHECK_CLOSE_REL(after.energy(), before.energy(), eps);
  // scattering has changed the direction
  BOOST_TEST(before.direction().dot(after.direction()) < 1);
  // scattering creates no new particles
  BOOST_TEST(outgoing.empty());
}
}  // namespace

BOOST_AUTO_TEST_SUITE(FatrasScattering)

BOOST_DATA_TEST_CASE(GeneralMixture, Dataset::parameters, pdg, phi, lambda, p,
                     seed) {
  test(GeneralMixtureScattering(), seed,
       Dataset::makeParticle(pdg, phi, lambda, p));
}

BOOST_DATA_TEST_CASE(GaussianMixture, Dataset::parameters, pdg, phi, lambda, p,
                     seed) {
  test(GaussianMixtureScattering(), seed,
       Dataset::makeParticle(pdg, phi, lambda, p));
}

BOOST_DATA_TEST_CASE(Highland, Dataset::parameters, pdg, phi, lambda, p, seed) {
  test(HighlandScattering(), seed, Dataset::makeParticle(pdg, phi, lambda, p));
}

BOOST_AUTO_TEST_SUITE_END()

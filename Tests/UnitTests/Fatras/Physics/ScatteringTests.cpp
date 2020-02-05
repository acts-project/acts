// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

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

constexpr double eps = 1e-10;

using GeneralMixtureScattering =
    ActsFatras::Scattering<ActsFatras::GeneralMixture>;
using GaussianMixtureScattering =
    ActsFatras::Scattering<ActsFatras::GaussianMixture>;
using HighlandScattering = ActsFatras::Scattering<ActsFatras::Highland>;

// Common test method that will be instantiated for each scattering model.
template <typename Scattering>
void run(const Scattering& scattering, const ActsFatras::Particle& before) {
  std::default_random_engine gen;
  ActsFatras::Particle after = before;

  const auto outgoing = scattering(gen, Dataset::thinSlab, after);
  // scattering leaves absolute energy/momentum unchanged
  CHECK_CLOSE_REL(after.momentum(), before.momentum(), eps);
  CHECK_CLOSE_REL(after.energy(), before.energy(), eps);
  // scattering creates no new particles
  BOOST_TEST(outgoing.empty());
}

}  // namespace

BOOST_AUTO_TEST_SUITE(FatrasScattering)

BOOST_DATA_TEST_CASE(GeneralMixture, Dataset::particleParameters, phi, lambda,
                     p, pdg, m, q) {
  run(GeneralMixtureScattering(),
      Dataset::makeParticle(phi, lambda, p, pdg, m, q));
}

BOOST_DATA_TEST_CASE(GaussianMixture, Dataset::particleParameters, phi, lambda,
                     p, pdg, m, q) {
  run(GaussianMixtureScattering(),
      Dataset::makeParticle(phi, lambda, p, pdg, m, q));
}

BOOST_DATA_TEST_CASE(Highland, Dataset::particleParameters, phi, lambda, p, pdg,
                     m, q) {
  run(HighlandScattering(), Dataset::makeParticle(phi, lambda, p, pdg, m, q));
}

BOOST_AUTO_TEST_SUITE_END()

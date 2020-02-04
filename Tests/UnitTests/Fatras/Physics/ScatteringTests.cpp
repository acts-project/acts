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
#include "ActsFatras/Physics/Scattering/Highland.hpp"
#include "ActsFatras/Physics/Scattering/Scattering.hpp"
#include "Dataset.hpp"

static constexpr double eps = 1e-10;

BOOST_AUTO_TEST_SUITE(FatrasScattering)

/// Test the scattering implementation
BOOST_DATA_TEST_CASE(HighlandScattering, Dataset::particleParameters, phi,
                     lambda, p, pdg, m, q) {
  using HighlandScattering = ActsFatras::Scattering<ActsFatras::Highland>;

  std::default_random_engine gen;
  ActsFatras::Particle before =
      Dataset::makeParticle(phi, lambda, p, pdg, m, q);
  ActsFatras::Particle after = before;

  HighlandScattering scattering;
  const auto outgoing = scattering(gen, Dataset::thinSlab, after);
  // scattering leaves absolute energy/momentum unchanged
  CHECK_CLOSE_REL(after.momentum(), before.momentum(), eps);
  CHECK_CLOSE_REL(after.energy(), before.energy(), eps);
  // scattering creates no new particles
  BOOST_TEST(outgoing.empty());
}

BOOST_AUTO_TEST_SUITE_END()

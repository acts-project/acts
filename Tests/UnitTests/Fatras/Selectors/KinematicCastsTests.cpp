// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "ActsFatras/Selectors/KinematicCasts.hpp"
#include "Dataset.hpp"

using namespace Acts::UnitLiterals;
using namespace ActsFatras::Casts;

static constexpr double eps = 1e-10;

BOOST_AUTO_TEST_SUITE(FatrasKinematicCasts)

BOOST_AUTO_TEST_CASE(CentralParticle) {
  const auto& particle = Dataset::centralPion;

  CHECK_SMALL(vR()(particle), eps);
  CHECK_SMALL(vZ()(particle), eps);
  CHECK_SMALL(AbsVz()(particle), eps);
  CHECK_SMALL(eta()(particle), eps);
  CHECK_SMALL(absEta()(particle), eps);
  CHECK_CLOSE_REL(pT()(particle), 1.5_GeV, eps);
  CHECK_CLOSE_REL(p()(particle), 1.5_GeV, eps);
  // TODO energy
}

BOOST_AUTO_TEST_CASE(ForwardParticle) {
  const auto& particle = Dataset::forwardPion;

  CHECK_SMALL(vR()(particle), eps);
  CHECK_CLOSE_REL(vZ()(particle), 100_mm, eps);
  CHECK_CLOSE_REL(AbsVz()(particle), 100_mm, eps);
  // TODO use value comparisons instead of relative checks
  BOOST_TEST(4.0 < eta()(particle));
  BOOST_TEST(4.0 < absEta()(particle));
  BOOST_TEST(10_MeV < pT()(particle));
  BOOST_TEST(1.5_GeV < p()(particle));
  // TODO energy
}

BOOST_AUTO_TEST_SUITE_END()

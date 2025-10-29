// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Seeding/CompositeSpacePointLineSeeder.hpp"

#include "StrawHitGeneratorHelper.hpp"

using namespace Acts;
using namespace Acts::Experimental;
using namespace Acts::Experimental::detail;
using namespace Acts::UnitLiterals;
using namespace Acts::detail::LineHelper;
using namespace Acts::PlanarHelper;
using namespace Acts::VectorHelpers;
using namespace Acts::Test;

using Seeder = CompositeSpacePointLineSeeder;

constexpr auto logLvl = Acts::Logging::Level::DEBUG;
constexpr std::size_t nEvents = 1;

ACTS_LOCAL_LOGGER(getDefaultLogger("StrawLineSeederTest", logLvl));

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(SeedingSuite)

BOOST_AUTO_TEST_CASE(CompositeSpacePointLineSeeder) {
  RandomEngine engine{1602};

  using GenCfg_t = MeasurementGenerator::Config;
  GenCfg_t genCfg{};
  genCfg.twinStraw = false;
  genCfg.createStrips = false;

  Seeder::Config seederCfg{};
  seederCfg.busyLayerLimit = 5;
  Seeder seeder{seederCfg};

  for (std::size_t evt = 0; evt < nEvents; ++evt) {
    ACTS_INFO("Generating event " << evt);
    const auto line = generateLine(engine, logger());
    auto testTubes =
        MeasurementGenerator::spawn(line, 0._ns, engine, genCfg, logger());
    std::unique_ptr<SpSorter> sorterPtr = std::make_unique<SpSorter>(testTubes);
    auto calibrator = std::make_unique<SpCalibrator>();

    using SeedOptions_t =
        Seeder::SeedOptions<Container_t, SpSorter, Container_t, SpCalibrator>;
    SeedOptions_t seedOpts{};
    seedOpts.splitter = std::move(sorterPtr);
    seedOpts.calibrator = calibrator.get();
    seedOpts.selector.connect<&isGoodHit>();
    ACTS_DEBUG(seedOpts);
    seeder.prepareSeedOptions(seedOpts);
    ACTS_DEBUG(seedOpts);

    while (auto seed = seeder.nextSeed(seedOpts)) {
      if (seed == std::nullopt)
        break;
      ACTS_INFO(*seed);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

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
#include "TFile.h"
#include "TTree.h"

using namespace Acts;
using namespace Acts::Experimental;
using namespace Acts::Experimental::detail;
using namespace Acts::UnitLiterals;
using namespace Acts::detail::LineHelper;
using namespace Acts::PlanarHelper;
using namespace Acts::VectorHelpers;
using namespace Acts::Test;

using Seeder = CompositeSpacePointLineSeeder;

constexpr auto logLvl = Acts::Logging::Level::INFO;
constexpr std::size_t nEvents = 5000;

ACTS_LOCAL_LOGGER(getDefaultLogger("StrawLineSeederTest", logLvl));

namespace ActsTests {

#define DECLARE_BRANCH(dTYPE, NAME) \
  dTYPE NAME{};                     \
  outTree->Branch(#NAME, &NAME);

BOOST_AUTO_TEST_SUITE(SeedingSuite)

void testSeeder(RandomEngine& engine, TFile& outFile) {
  auto outTree = std::make_unique<TTree>("StrawSeederTree", "StrawSeederTree");

  DECLARE_BRANCH(double, trueY0);
  DECLARE_BRANCH(double, trueTheta);

  DECLARE_BRANCH(std::vector<double>, recoY0);
  DECLARE_BRANCH(std::vector<double>, recoTheta);
  DECLARE_BRANCH(std::vector<double>, uncertY0);
  DECLARE_BRANCH(std::vector<double>, uncertTheta);
  DECLARE_BRANCH(std::vector<uint>, nStraws);
  DECLARE_BRANCH(std::vector<uint>, nStrips);
  DECLARE_BRANCH(uint, nSeeds);

  using GenCfg_t = MeasurementGenerator::Config;
  GenCfg_t genCfg{};
  genCfg.twinStraw = false;
  genCfg.createStrips = false;

  Seeder::Config seederCfg{};
  seederCfg.busyLayerLimit = 20;
  Seeder seeder{seederCfg};

  for (std::size_t evt = 0; evt < nEvents; ++evt) {
    if(evt%100 == 0) ACTS_INFO("Generating event " << evt);
    const auto line = generateLine(engine, logger());
    auto linePars = line.parameters();
    trueY0 = linePars[toUnderlying(Line_t::ParIndex::y0)];
    trueTheta = linePars[toUnderlying(Line_t::ParIndex::theta)];
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
    ACTS_DEBUG("Preparing seed options ");
    seeder.prepareSeedOptions(seedOpts);
    ACTS_DEBUG(seedOpts);
    nSeeds = 0;
    while (auto seed = seeder.nextSeed(seedOpts)) {
      ACTS_DEBUG("Seed finder loop "<< seedOpts);
      if (seed == std::nullopt)
        break;
      recoY0.push_back(seed->y0);
      recoTheta.push_back(seed->theta);
      uncertY0.push_back(seed->dY0);
      uncertTheta.push_back(seed->dTheta);
      nStraws.push_back(seed->nStrawHits);
      nSeeds++;
    }

    outTree->Fill();
    recoY0.clear();
    recoTheta.clear();
    uncertY0.clear();
    uncertTheta.clear();
    nStraws.clear();
    nStrips.clear();
  }
  outFile.WriteObject(outTree.get(), outTree->GetName());
}

BOOST_AUTO_TEST_CASE(CompositeSpacePointLineSeeder) {
  RandomEngine engine{1602};
  std::unique_ptr<TFile> outFile{
      TFile::Open("StrawLineSeedTest.root", "RECREATE")};

  BOOST_CHECK_EQUAL(outFile->IsZombie(), false);

  testSeeder(engine, *outFile);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

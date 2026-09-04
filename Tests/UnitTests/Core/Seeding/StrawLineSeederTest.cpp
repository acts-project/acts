// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <set>

#include "StrawHitGeneratorHelper.hpp"
#include "TFile.h"
#include "TTree.h"

constexpr auto logLvl = Acts::Logging::Level::INFO;
constexpr std::size_t nEvents = 5;

ACTS_LOCAL_LOGGER(getDefaultLogger("StrawLineSeederTest", logLvl));

#define DECLARE_BRANCH(dTYPE, NAME) \
  dTYPE NAME{};                     \
  outTree->Branch(#NAME, &NAME);

namespace ActsTests {

void testSeeder(RandomEngine& engine, TFile& outFile) {
  auto outTree = std::make_unique<TTree>("StrawSeederTree", "StrawSeederTree");

  DECLARE_BRANCH(double, trueY0);
  DECLARE_BRANCH(double, trueTheta);
  DECLARE_BRANCH(std::size_t, nTruthStraws);
  DECLARE_BRANCH(std::size_t, nTruthStrips);
  DECLARE_BRANCH(std::vector<double>, recoY0);
  DECLARE_BRANCH(std::vector<double>, recoTheta);
  DECLARE_BRANCH(std::vector<std::size_t>, nStraws);
  DECLARE_BRANCH(std::vector<std::size_t>, nStrips);
  DECLARE_BRANCH(std::size_t, nSeeds);

  using GenCfg_t = MeasurementGenerator::Config;
  GenCfg_t genCfg{};
  genCfg.twinStraw = false;
  genCfg.createStrips = true;

  CompositeSpacePointLineSeeder::Config seederCfg{};
  seederCfg.busyLayerLimit = 20;
  const CompositeSpacePointLineSeeder seeder{
      seederCfg, getDefaultLogger("StrawLineSeederTest", logLvl)};

  for (std::size_t evt = 0; evt < nEvents; ++evt) {
    if (evt % 100 == 0) {
      ACTS_INFO("Generating event " << evt);
    }
    const auto line = generateLine(engine, logger());
    auto linePars = line.parameters();
    trueY0 = linePars[toUnderlying(FitParIndex::y0)];
    trueTheta = linePars[toUnderlying(FitParIndex::theta)];
    auto testTubes =
        MeasurementGenerator::spawn(line, 0._ns, engine, genCfg, logger());
    nTruthStraws = static_cast<std::size_t>(std::ranges::count_if(
        testTubes, [](const auto& testSp) { return testSp->isStraw(); }));
    nTruthStrips = static_cast<std::size_t>(std::ranges::count_if(
        testTubes, [](const auto& testSp) { return !testSp->isStraw(); }));
    auto calibrator = std::make_unique<SpCalibrator>();

    using SeedState_t =
        CompositeSpacePointLineSeeder::SeedingState<Container_t, Container_t,
                                                    SpSorter>;
    SeedState_t seedState{startParameters(line, testTubes), testTubes,
                          calibrator.get()};
    ACTS_DEBUG(seedState);
    nSeeds = 0;
    CalibrationContext cctx{};
    while (auto seed = seeder.nextSeed(cctx, seedState)) {
      ACTS_DEBUG("Seed finder loop " << seedState);
      if (seed == std::nullopt) {
        break;
      }
      recoTheta.push_back(seed->parameters[toUnderlying(FitParIndex::theta)]);
      recoY0.push_back(seed->parameters[toUnderlying(FitParIndex::y0)]);

      const auto seedStraws = static_cast<std::size_t>(std::ranges::count_if(
          seed->hits, [](const auto& testMe) { return testMe->isStraw(); }));
      const auto seedStrips = seed->hits.size() - seedStraws;
      nStraws.push_back(seedStraws);
      nStrips.push_back(seedStrips);
      ++nSeeds;
    }
    ACTS_DEBUG("======Event " << evt << " found " << nSeeds << " seeds.");

    outTree->Fill();
    recoY0.clear();
    recoTheta.clear();
    nStraws.clear();
    nStrips.clear();
  }
  outFile.WriteObject(outTree.get(), outTree->GetName());
}

using GenCfg_t = MeasurementGenerator::Config;

BOOST_AUTO_TEST_SUITE(SeedingSuite)
BOOST_AUTO_TEST_CASE(SeederTest) {
  RandomEngine engine{1602};
  std::unique_ptr<TFile> outFile{
      TFile::Open("StrawLineSeedTest.root", "RECREATE")};

  BOOST_CHECK_EQUAL(outFile->IsZombie(), false);

  testSeeder(engine, *outFile);
}
BOOST_AUTO_TEST_CASE(SeedTangents) {
  RandomEngine engine{117};
  constexpr double tolerance = 1.e-3;

  using Seeder = CompositeSpacePointLineSeeder;
  using SeedAux = CompSpacePointAuxiliaries;
  using enum Seeder::TangentAmbi;
  GenCfg_t genCfg{};
  genCfg.createStrips = false;
  genCfg.createStraws = true;
  genCfg.smearRadius = false;
  for (std::size_t evt = 0; evt < nEvents; ++evt) {
    const auto line = generateLine(engine, logger());
    auto testTubes =
        MeasurementGenerator::spawn(line, 0._ns, engine, genCfg, logger());

    const double lineTanBeta = line.direction().y() / line.direction().z();
    const double lineY0 = line.position().y();

    for (std::size_t m1 = testTubes.size() - 1; m1 > testTubes.size() / 2;
         --m1) {
      for (std::size_t m2 = 0; m2 < m1; ++m2) {
        const auto& bottomTube = *testTubes[m2];
        const auto& topTube = *testTubes[m1];

        const int signTop = SeedAux::strawSign(line, topTube);
        const int signBot = SeedAux::strawSign(line, bottomTube);
        const auto trueAmbi = Seeder::encodeAmbiguity(signTop, signBot);

        ACTS_DEBUG(__func__ << "() " << __LINE__ << " - "
                            << std::format("bottom tube @ {:}, r: {:.3f}({:})",
                                           toString(bottomTube.localPosition()),
                                           (signBot * bottomTube.driftRadius()),
                                           signBot > 0 ? "R" : "L")
                            << ", "
                            << std::format("top tube @ {:}, r: {:.3f} ({:})",
                                           toString(topTube.localPosition()),
                                           (signTop * topTube.driftRadius()),
                                           signTop > 0 ? "R" : "L"));

        bool seenTruePars{false};
        std::set<std::pair<double, double>> fourSeedPars{};
        for (const auto ambi : {LL, RL, LR, RR}) {
          const auto pars =
              Seeder::constructTangentLine(topTube, bottomTube, ambi);

          const bool isTruePars =
              Acts::abs(std::tan(pars.theta) - lineTanBeta) < tolerance &&
              Acts::abs(lineY0 - pars.y0) < tolerance;
          seenTruePars |= isTruePars;
          ACTS_VERBOSE(__func__
                       << "() " << __LINE__ << " - Test ambiguity "
                       << CompositeSpacePointLineSeeder::toString(ambi)
                       << " -> "
                       << std::format("theta: {:.3f}, ", pars.theta / 1._degree)
                       << std::format("tanTheta: {:.3f}, ",
                                      std::tan(pars.theta))
                       << std::format("y0: {:.3f}", pars.y0)
                       << (!isTruePars ? (ambi == trueAmbi ? " xxxxxx" : "")
                                       : " <-------"));
          const Vector3 seedPos = pars.y0 * Vector3::UnitY();
          const Vector3 seedDir =
              makeDirectionFromAxisTangents(0., std::tan(pars.theta));

          const double chi2Top = SeedAux::chi2Term(seedPos, seedDir, topTube);
          const double chi2Bot =
              SeedAux::chi2Term(seedPos, seedDir, bottomTube);
          ACTS_VERBOSE(__func__ << "() " << __LINE__
                                << " - Resulting chi2 top: " << chi2Top
                                << ", bottom: " << chi2Bot);
          BOOST_CHECK_LE(chi2Top, 1.e-17);
          BOOST_CHECK_LE(chi2Bot, 1.e-17);
          BOOST_CHECK_EQUAL(
              fourSeedPars.insert(std::make_pair(pars.theta, pars.y0)).second,
              true);
        }
        BOOST_CHECK_EQUAL(seenTruePars, true);
        BOOST_CHECK_EQUAL(fourSeedPars.size(), 4);
        const auto seedPars =
            Seeder::constructTangentLine(topTube, bottomTube, trueAmbi);
        /// Construct line parameters
        ACTS_DEBUG(__func__
                   << "() " << __LINE__ << " - Line tan beta: " << lineTanBeta
                   << ", reconstructed tan theta: " << std::tan(seedPars.theta)
                   << ", line y0: " << lineY0
                   << ", reconstructed y0: " << seedPars.y0);
        BOOST_CHECK_CLOSE(std::tan(seedPars.theta), lineTanBeta, tolerance);
        BOOST_CHECK_CLOSE(seedPars.y0, lineY0, tolerance);
      }
    }
  }
}
}  // namespace ActsTests
BOOST_AUTO_TEST_SUITE_END()

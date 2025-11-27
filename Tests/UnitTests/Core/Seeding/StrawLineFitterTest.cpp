// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <chrono>
#include <format>
#include <random>
#include <ranges>
#include <set>
#include <span>

#include <TFile.h>
#include <TTree.h>

#include "StrawHitGeneratorHelper.hpp"

using namespace Acts;
using namespace Acts::Experimental;
using namespace Acts::Experimental::detail;
using namespace Acts::UnitLiterals;
using namespace Acts::detail::LineHelper;
using namespace Acts::PlanarHelper;
using namespace Acts::VectorHelpers;

using TimePoint_t = std::chrono::system_clock::time_point;
using RandomEngine = std::mt19937;
using uniform = std::uniform_real_distribution<double>;
using normal_t = std::normal_distribution<double>;
using Line_t = CompSpacePointAuxiliaries::Line_t;
using ResidualIdx = CompSpacePointAuxiliaries::ResidualIdx;
using FitParIndex = CompSpacePointAuxiliaries::FitParIndex;
using ParamVec_t = CompositeSpacePointLineFitter::ParamVec_t;
using Fitter_t = CompositeSpacePointLineFitter;

constexpr auto logLvl = Acts::Logging::Level::INFO;
constexpr std::size_t nEvents = 200000;

ACTS_LOCAL_LOGGER(getDefaultLogger("StrawLineFitterTest", logLvl));

namespace ActsTests {

using GenCfg_t = MeasurementGenerator::Config;

BOOST_AUTO_TEST_SUITE(SeedingSuite)

BOOST_AUTO_TEST_CASE(SeedTangents) {
  RandomEngine engine{1602};
  constexpr double tolerance = 1.e-3;

  using Seeder = CompositeSpacePointLineSeeder;
  using SeedAux = CompSpacePointAuxiliaries;
  using enum Seeder::TangentAmbi;
  GenCfg_t genCfg{};
  genCfg.createStrips = false;
  genCfg.createStraws = true;
  genCfg.smearRadius = false;

  for (std::size_t evt = 0; evt < 100; ++evt) {
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
                   << "() " << __LINE__ << " - Line tan theta: " << lineTanBeta
                   << ", reconstructed theta: " << std::tan(seedPars.theta)
                   << ", line y0: " << lineY0
                   << ", reconstructed y0: " << seedPars.y0);
        BOOST_CHECK_CLOSE(std::tan(seedPars.theta), lineTanBeta, tolerance);
        BOOST_CHECK_CLOSE(seedPars.y0, lineY0, tolerance);
      }
    }
  }
}

#define DECLARE_BRANCH(dType, bName) \
  dType bName{};                     \
  outTree->Branch(#bName, &bName);

void runFitTest(const Fitter_t::Config& fitCfg, const GenCfg_t& genCfg,
                const std::string& testName, RandomEngine& engine,
                TFile& outFile) {
  Fitter_t fitter{fitCfg, getDefaultLogger(
                              std::format("LineFitter_{:}", testName), logLvl)};

  auto outTree = std::make_unique<TTree>(
      std::format("{:}Tree", testName).c_str(), "MonitorTree");

  TimePoint_t start = std::chrono::system_clock::now();

  ACTS_INFO("Start test " << testName << ".");

  DECLARE_BRANCH(double, trueY0);
  DECLARE_BRANCH(double, trueX0);
  DECLARE_BRANCH(double, trueTheta);
  DECLARE_BRANCH(double, truePhi);
  DECLARE_BRANCH(double, trueT0);
  DECLARE_BRANCH(double, trueProjTheta);
  DECLARE_BRANCH(double, trueProjPhi);

  DECLARE_BRANCH(double, recoY0);
  DECLARE_BRANCH(double, recoX0);
  DECLARE_BRANCH(double, recoTheta);
  DECLARE_BRANCH(double, recoPhi);
  DECLARE_BRANCH(double, recoT0);
  DECLARE_BRANCH(double, recoProjTheta);
  DECLARE_BRANCH(double, recoProjPhi);

  DECLARE_BRANCH(double, sigmaY0);
  DECLARE_BRANCH(double, sigmaX0);
  DECLARE_BRANCH(double, sigmaTheta);
  DECLARE_BRANCH(double, sigmaPhi);
  DECLARE_BRANCH(double, sigmaT0);

  DECLARE_BRANCH(double, chi2);
  DECLARE_BRANCH(unsigned, nIter);
  DECLARE_BRANCH(unsigned, nDoF);

  /// @brief Fill the parameter array to the tree variables
  /// @param pars: Parameter array to safe
  /// @param y0: Reference to the variable storing y0
  /// @param x0: Reference to the variable storing x0
  /// @param theta: Reference to the variable storing theta
  /// @param phi: Reference to the variable storing phi
  auto fillPars = [](const auto pars, double& y0, double& x0, double& theta,
                     double& phi) {
    y0 = pars[toUnderlying(FitParIndex::y0)];
    x0 = pars[toUnderlying(FitParIndex::x0)];
    theta = pars[toUnderlying(FitParIndex::theta)] / 1._degree;
    phi = pars[toUnderlying(FitParIndex::phi)] / 1._degree;
  };
  /// @brief Fill the
  auto fillProjected = [](const auto pars, double& projTheta, double& projPhi) {
    auto dir =
        makeDirectionFromPhiTheta(pars[toUnderlying(FitParIndex::phi)],
                                  pars[toUnderlying(FitParIndex::theta)]);
    projTheta = std::atan(dir[ePos1] / dir[ePos2]) / 1._degree;
    projPhi = std::atan(dir[ePos0] / dir[ePos2]) / 1._degree;
  };
  // Pass a localToGlobal transform to the calibrator to proper handling the ToF
  auto calibrator = std::make_unique<SpCalibrator>();
  for (std::size_t evt = 0; evt < nEvents; ++evt) {
    const auto line = generateLine(engine, logger());
    fillPars(line.parameters(), trueY0, trueX0, trueTheta, truePhi);
    fillProjected(line.parameters(), trueProjTheta, trueProjPhi);
    const double t0 = uniform{-50._ns, 50._ns}(engine);
    trueT0 = t0 / 1._ns;

    using FitOpts_t = Fitter_t::FitOptions<Container_t, SpCalibrator>;

    FitOpts_t fitOpts{};
    fitOpts.calibrator = calibrator.get();

    fitOpts.selector.connect<&isGoodHit>();
    fitOpts.measurements =
        MeasurementGenerator::spawn(line, t0, engine, genCfg, logger());
    fitOpts.startParameters = startParameters(line, fitOpts.measurements);
    fillPars(fitOpts.startParameters, recoY0, recoX0, recoTheta, recoPhi);

    //

    auto result = fitter.fit(std::move(fitOpts));
    if (!result.converged) {
      ACTS_INFO("Fit failed - not converged.");
      continue;
    }
    ACTS_DEBUG("Fit Successful.");

    fillPars(result.parameters, recoY0, recoX0, recoTheta, recoPhi);
    fillProjected(result.parameters, recoProjTheta, recoProjPhi);

    recoT0 = result.parameters[toUnderlying(FitParIndex::t0)] / 1._ns;

    auto extractUncert = [&result](const auto idx) {
      return std::sqrt(result.covariance(toUnderlying(idx), toUnderlying(idx)));
    };
    sigmaY0 = extractUncert(FitParIndex::y0);
    sigmaX0 = extractUncert(FitParIndex::x0);
    sigmaTheta = extractUncert(FitParIndex::theta) / 1._degree;
    sigmaPhi = extractUncert(FitParIndex::phi) / 1._degree;
    sigmaT0 = extractUncert(FitParIndex::t0) / 1._ns;

    chi2 = result.chi2;
    nDoF = result.nDoF;
    nIter = result.nIter;

    outTree->Fill();
    if ((evt + 1) % 10 == 0u) {
      ACTS_INFO("Processed " << (evt + 1) << "/" << nEvents << " events.");
    }
  }

  TimePoint_t end = std::chrono::system_clock::now();  // timing: get end time
  auto diff =
      std::chrono::duration_cast<std::chrono::seconds>(end - start).count();

  outFile.WriteObject(outTree.get(), outTree->GetName());
  ACTS_INFO("Test finished. " << outTree->GetEntries()
                              << " tracks written. Test took " << diff
                              << " seconds.");
}
#undef DECLARE_BRANCH

BOOST_AUTO_TEST_CASE(SimpleLineFit) {
  auto outFile =
      std::make_unique<TFile>("StrawLineFitterTest.root", "RECREATE");

  Fitter_t::Config fitCfg{};
  fitCfg.useHessian = true;
  fitCfg.calcAlongStraw = true;
  fitCfg.recalibrate = true;

  GenCfg_t genCfg{};
  genCfg.twinStraw = false;
  genCfg.createStrips = false;
  // 2D straw only test
  {
    RandomEngine engine{1602};
    runFitTest(fitCfg, genCfg, "StrawOnlyTest", engine, *outFile);
  }
  // 2D straw only test with time
  {
    fitCfg.fitT0 = true;
    fitCfg.useFastFitter = true;
    fitCfg.fastPreFitter = true;
    fitCfg.precCutOff = 1e-4;
    RandomEngine engine{1602};
    runFitTest(fitCfg, genCfg, "StrawOnlyTestT0", engine, *outFile);
  }
  fitCfg.fitT0 = false;
  // fast straw only test
  {
    fitCfg.useFastFitter = true;
    fitCfg.fastPreFitter = false;
    RandomEngine engine{1503};
    runFitTest(fitCfg, genCfg, "FastStrawOnlyTest", engine, *outFile);
  }
  // 2D straws + twin measurement test
  {
    fitCfg.useFastFitter = false;
    genCfg.twinStraw = true;
    RandomEngine engine{1701};
    runFitTest(fitCfg, genCfg, "StrawAndTwinTest", engine, *outFile);
  }
  genCfg.createStrips = true;
  genCfg.twinStraw = false;
  genCfg.combineSpacePoints = false;
  // 1D straws + single strip measurements
  {
    RandomEngine engine{1404};
    runFitTest(fitCfg, genCfg, "StrawAndStripTest", engine, *outFile);
  }
  // Strip only
  {
    genCfg.createStraws = false;
    genCfg.combineSpacePoints = true;
    genCfg.stripPitchLoc1 = 500._um;
    RandomEngine engine{2070};
    runFitTest(fitCfg, genCfg, "StripOnlyTest", engine, *outFile);
  }
  // Strip stereo test
  {
    genCfg.stripDirLoc1 = makeDirectionFromPhiTheta(60._degree, 90._degree);
    RandomEngine engine{2225};
    fitCfg.parsToUse = {FitParIndex::x0, FitParIndex::y0, FitParIndex::theta,
                        FitParIndex::phi};
    runFitTest(fitCfg, genCfg, "StereoStripTest", engine, *outFile);
  }
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace ActsTests

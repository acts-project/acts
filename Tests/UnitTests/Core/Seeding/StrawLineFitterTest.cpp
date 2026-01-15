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
#include <future>
#include <iostream>
#include <random>
#include <ranges>
#include <span>
#include <thread>

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

#include "StrawHitGeneratorHelper.hpp"

using TimePoint_t = std::chrono::system_clock::time_point;
using Fitter_t = CompositeSpacePointLineFitter;

constexpr auto logLvl = Acts::Logging::Level::INFO;
constexpr std::size_t nEvents = 1;
constexpr long int nThreads = 1;
std::mutex writeMutex{};

ACTS_LOCAL_LOGGER(getDefaultLogger("StrawLineFitterTest", logLvl));

namespace ActsTests {

using GenCfg_t = MeasurementGenerator::Config;

#define DECLARE_BRANCH(dType, bName) \
  dType bName{};                     \
  outTree->Branch(#bName, &bName);
// NOLINTBEGIN
long int runFitTest(Fitter_t::Config fitCfg, GenCfg_t genCfg,
                    const std::string& testName, const unsigned seed,
                    TFile& outFile) {
  // NOLINTEND
  auto outTree = std::make_unique<TTree>(
      std::format("{:}Tree", testName).c_str(), "MonitorTree");
  outTree->SetDirectory(nullptr);

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
  DECLARE_BRANCH(char, converged);

  RandomEngine engine{seed};
  Fitter_t fitter{fitCfg, getDefaultLogger(
                              std::format("LineFitter_{:}", testName), logLvl)};

  TimePoint_t start = std::chrono::system_clock::now();

  ACTS_INFO("Start test " << testName << ".");

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
    Vector3 properDir =
        makeDirectionFromPhiTheta(pars[toUnderlying(FitParIndex::phi)],
                                  pars[toUnderlying(FitParIndex::theta)]);
    theta =
        VectorHelpers::theta(copySign(properDir, properDir.z())) / 1._degree;
    phi = VectorHelpers::phi(copySign(properDir, properDir.z())) / 1._degree;
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
  std::size_t goodFits{0};
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
      ACTS_DEBUG("Fit " << outTree->GetName() << " failed.");
      converged = 0;
      chi2 = -1.;
      nDoF = 1;
      nIter = fitter.config().maxIter;
      outTree->Fill();
      continue;
    }

    ACTS_DEBUG("Fit Successful.");
    converged = 1;
    ++goodFits;
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
    if ((evt + 1) % 10000 == 0u) {
      ACTS_INFO("Processed " << (evt + 1) << "/" << nEvents
                             << " events. Test: " << outTree->GetName());
    }
  }

  TimePoint_t end = std::chrono::system_clock::now();  // timing: get end time
  auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                  .count();
  std::unique_lock guard{writeMutex};
  outFile.WriteObject(outTree.get(), outTree->GetName());
  ACTS_INFO("Test " << outTree->GetName() << " finished. " << goodFits
                    << " tracks written. It took " << (diff / 1000)
                    << " seconds.");
  return diff;
}
#undef DECLARE_BRANCH

BOOST_AUTO_TEST_SUITE(StrawLineFitTestSuite)

BOOST_AUTO_TEST_CASE(SimpleLineFit) {
  using namespace std::chrono_literals;

  auto outFile =
      std::make_unique<TFile>("StrawLineFitterTest.root", "RECREATE");

  Fitter_t::Config fitCfg{};
  fitCfg.useHessian = false;
  fitCfg.calcAlongStraw = true;
  fitCfg.recalibrate = false;
  fitCfg.useFastFitter = false;
  fitCfg.ranges[toUnderlying(FitParIndex::theta)] =
      std::array{1._degree, 179._degree};
  fitCfg.ranges[toUnderlying(FitParIndex::phi)] =
      std::array{-179._degree, 179._degree};
  fitCfg.ranges[toUnderlying(FitParIndex::x0)] = std::array{-1000., 1000.};
  fitCfg.ranges[toUnderlying(FitParIndex::y0)] = std::array{-1000., 1000.};
  /// Configuration for fast pre-fit
  Fitter_t::Config fastPreCfg{fitCfg};
  fastPreCfg.useFastFitter = true;
  /// Configuration for fast only fit
  Fitter_t::Config fastCfg{fastPreCfg};
  fastCfg.fastPreFitter = false;

  // 2D straw only test
  std::vector<std::pair<std::string, std::future<long int>>> timings{};
  using namespace std::chrono_literals;

  auto sendSleep = [&timings]() {
    do {
      std::this_thread::sleep_for(100ms);
    } while (std::ranges::count_if(timings, [](const auto& timeObj) {
               return timeObj.second.wait_for(0ms) != std::future_status::ready;
             }) >= nThreads);
  };
  /// @brief Helper lambda to launch the test for a given measurement configuration.
  ///        Always three tests are launched:
  ///           Fast: Use only the fast fitter to obtain the results
  ///           FastPre: Use the fast fitter for a pre estimate of the
  ///           parameters
  ///                    followed by the full fitter
  ///           Full: Just use the full fitter to obtain the results
  /// @param testName: Name of the tests
  /// @param genCfg: Configuration object to generate the measurements
  /// @param seed: Seed number for the random number generator
  auto launchTest = [&](const std::string& testName, const GenCfg_t& genCfg,
                        const unsigned seed) {
    sendSleep();
    timings.emplace_back(
        "Fast" + testName, std::async(std::launch::async, [&]() {
          return runFitTest(fastCfg, genCfg, "Fast" + testName, seed, *outFile);
        }));
    sendSleep();
    timings.emplace_back(testName, std::async(std::launch::async, [&]() {
                           return runFitTest(fitCfg, genCfg, testName, seed,
                                             *outFile);
                         }));

    sendSleep();
    timings.emplace_back(
        "FastPre" + testName, std::async(std::launch::async, [&]() {
          return runFitTest(fastPreCfg, genCfg, "FastPre" + testName, seed,
                            *outFile);
        }));
    sendSleep();
  };
  {
    GenCfg_t genCfg{};
    genCfg.twinStraw = false;
    genCfg.createStrips = false;
    launchTest("StrawOnlyTest", genCfg, 1602);
  }
  // 2D straws + twin measurement test
  {
    GenCfg_t genCfg{};
    genCfg.createStraws = true;
    genCfg.twinStraw = true;
    genCfg.createStrips = false;

    launchTest("StrawAndTwinTest", genCfg, 1503);
  }
  // 1D straws + single strip measurements
  {
    GenCfg_t genCfg{};
    genCfg.createStrips = true;
    genCfg.twinStraw = false;
    genCfg.combineSpacePoints = false;
    genCfg.discretizeStrips = true;
    genCfg.createStraws = true;
    genCfg.stripPitchLoc1 = 2._cm;
    genCfg.stripPitchLoc0 = 3.5_cm;
    launchTest("StrawAndStripTest", genCfg, 1701);
  }
  // 1D straws + 2D strip measurements
  {
    GenCfg_t genCfg{};

    genCfg.createStrips = true;
    genCfg.twinStraw = false;
    genCfg.combineSpacePoints = true;
    genCfg.discretizeStrips = true;
    genCfg.createStraws = true;
    genCfg.stripPitchLoc1 = 2._cm;
    genCfg.stripPitchLoc0 = 3.5_cm;
    launchTest("StrawAndStrip2DTest", genCfg, 1404);
  }
  // Strip only
  {
    GenCfg_t genCfg{};
    genCfg.createStrips = true;
    genCfg.twinStraw = false;
    genCfg.combineSpacePoints = false;
    genCfg.discretizeStrips = true;
    genCfg.createStraws = false;
    genCfg.stripPitchLoc1 = 500._um;
    genCfg.stripPitchLoc0 = 3._cm;
    launchTest("StripOnlyTest", genCfg, 2070);
  }
  // 2D Strip only
  {
    GenCfg_t genCfg{};
    genCfg.createStrips = true;
    genCfg.twinStraw = false;
    genCfg.combineSpacePoints = true;
    genCfg.discretizeStrips = true;
    genCfg.createStraws = false;
    genCfg.stripPitchLoc1 = 500._um;
    genCfg.stripPitchLoc0 = 3._cm;
    launchTest("Strip2DOnlyTest", genCfg, 2225);
  }
  // Strip stereo test
  {
    GenCfg_t genCfg{};
    genCfg.createStrips = true;
    genCfg.twinStraw = false;
    genCfg.combineSpacePoints = false;
    genCfg.discretizeStrips = true;
    genCfg.createStraws = false;
    genCfg.stripPitchLoc1 = 500._um;
    genCfg.stripDirLoc0.clear();
    genCfg.stripDirLoc1 = {
        makeDirectionFromPhiTheta(0._degree, 90._degree),
        makeDirectionFromPhiTheta(0._degree, 90._degree),
        makeDirectionFromPhiTheta(-1.5_degree, 90._degree),
        makeDirectionFromPhiTheta(1.5_degree, 90._degree),
        makeDirectionFromPhiTheta(-1.5_degree, 90._degree),
        makeDirectionFromPhiTheta(1.5_degree, 90._degree),
        makeDirectionFromPhiTheta(0._degree, 90._degree),
        makeDirectionFromPhiTheta(0._degree, 90._degree),

    };

    fitCfg.parsToUse = {FitParIndex::x0, FitParIndex::y0, FitParIndex::theta,
                        FitParIndex::phi};
    launchTest("StereoStripTest", genCfg, 1800);
  }
  /// Wait until all tests ar ecompleted
  while (std::ranges::any_of(timings, [](const auto& lblTh) {
    return lblTh.second.wait_for(0ms) != std::future_status::ready;
  })) {
    std::this_thread::sleep_for(10ms);
  }
  {
    auto timeHisto = std::make_unique<TH1D>("TestTimings", "timings",
                                            timings.size(), 0, timings.size());
    int bin{1};
    for (auto& [label, result] : timings) {
      timeHisto->GetXaxis()->SetBinLabel(bin, label.c_str());
      timeHisto->SetBinContent(bin, static_cast<double>(result.get()) /
                                        static_cast<double>(nEvents));
      ++bin;
    }
    outFile->WriteObject(timeHisto.get(), timeHisto->GetName());
  }
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace ActsTests

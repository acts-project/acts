// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsExamples/Validation/PatternRecognitionPerformanceCollector.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <numbers>
#include <random>
#include <vector>

using namespace Acts;
using namespace ActsExamples;

namespace {

void checkAgainstBruteForce(const std::vector<std::array<double, 2>>& points,
                            double maxDeltaR = 0.3) {
  SimParticleContainer particles;
  ParticleTrackMatching matching;
  TrackParticleMatching trackMatching;
  auto trackBackend = std::make_shared<VectorTrackContainer>();
  auto stateBackend = std::make_shared<VectorMultiTrajectory>();
  ActsExamples::TrackContainer tracks(trackBackend, stateBackend);
  for (std::size_t i = 0; i < points.size(); ++i) {
    auto barcode = SimBarcode{}.withVertexPrimary(1).withParticle(i + 1);
    SimParticle particle(barcode, PdgParticle::eMuon);
    const auto [eta, phi] = points[i];
    particle.initialState().setDirection(
        std::isinf(eta)
            ? Vector3(0., 0., std::copysign(1., eta))
            : Vector3(std::cos(phi), std::sin(phi), std::sinh(eta)));
    particle.initialState().setAbsoluteMomentum(10.);
    particles.insert(particle);
    if (i % 2 == 0) {
      auto track = tracks.makeTrack();
      matching[barcode].track = track.index();
      trackMatching[track.index()] = {
          TrackMatchClassification::Matched, barcode, {}};
    }
  }
  ConstTrackContainer constTracks(
      std::make_shared<ConstVectorTrackContainer>(*trackBackend),
      std::make_shared<ConstVectorMultiTrajectory>(*stateBackend));
  PatternRecognitionPerformanceCollector::Config config;
  config.effPlotToolConfig.varBinning["DeltaR"] =
      EffPlotTool::BoostRegularAxis(100, 0., maxDeltaR, "deltaR");
  PatternRecognitionPerformanceCollector collector(
      config, getDefaultLogger("Test", Logging::WARNING));
  EffPlotTool reference(config.effPlotToolConfig, Logging::WARNING);
  auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  collector.fill(gctx, constTracks, particles, trackMatching, matching, {});
  for (const auto& particle : particles) {
    double nearest = -1.;
    for (const auto& other : particles) {
      if (particle.particleId() == other.particleId()) {
        continue;
      }
      const double distance =
          VectorHelpers::deltaR(particle.direction(), other.direction());
      if (nearest == -1. || distance < nearest) {
        nearest = distance;
      }
    }
    reference.fill(gctx, particle.initialState(), nearest,
                   matching.contains(particle.particleId()));
  }
  for (const auto& [name, expected] : reference.efficiencies1D()) {
    BOOST_TEST_CONTEXT(name) {
      const auto& actual = collector.effPlotTool().efficiencies1D().at(name);
      for (int bin = -1; bin <= expected.totalHistogram().axis(0).size();
           ++bin) {
        BOOST_CHECK_EQUAL(actual.totalHistogram().at(bin).value(),
                          expected.totalHistogram().at(bin).value());
        BOOST_CHECK_EQUAL(actual.acceptedHistogram().at(bin).value(),
                          expected.acceptedHistogram().at(bin).value());
      }
    }
  }
}

}  // namespace

BOOST_AUTO_TEST_CASE(EmptySingletonAndIsolatedParticles) {
  checkAgainstBruteForce({});
  checkAgainstBruteForce({{0., 0.}});
  checkAgainstBruteForce({{-5., -2.}, {0., 0.}, {5., 2.}});
}

BOOST_AUTO_TEST_CASE(PeriodicPhiCoincidenceAndCutoff) {
  const double pi = std::numbers::pi;
  checkAgainstBruteForce({{0., pi - 0.001},
                          {0., -pi + 0.001},
                          {1., 0.},
                          {1., 0.},
                          {2., 0.},
                          {2., 0.3},
                          {3., 0.},
                          {3., std::nextafter(0.3, 0.)},
                          {4., 0.},
                          {4., std::nextafter(0.3, 1.)}});
}

BOOST_AUTO_TEST_CASE(ConfiguredRangeAndRandomParticles) {
  std::mt19937 rng(42);
  std::uniform_real_distribution<double> eta(-6., 6.);
  std::uniform_real_distribution<double> phi(-std::numbers::pi,
                                             std::numbers::pi);
  std::vector<std::array<double, 2>> points;
  for (std::size_t i = 0; i < 512; ++i) {
    points.push_back({eta(rng), phi(rng)});
  }
  for (double maxDeltaR : {0.01, 0.3, 1., 10.}) {
    BOOST_TEST_CONTEXT(maxDeltaR) {
      checkAgainstBruteForce(points, maxDeltaR);
    }
  }
}

BOOST_AUTO_TEST_CASE(BeamAxisParticles) {
  const double infinity = std::numeric_limits<double>::infinity();
  checkAgainstBruteForce({{infinity, 0.}, {0., 0.}, {0., 0.1}});
  checkAgainstBruteForce({{-infinity, 0.}, {0., 0.}});
}

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AnnealingUtility.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/AdaptiveGridDensityVertexFinder.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFinder.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFitter.hpp"
#include "Acts/Vertexing/GaussianTrackDensity.hpp"
#include "Acts/Vertexing/GridDensityVertexFinder.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/IVertexFinder.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/TrackDensityVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingError.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>
#include <numbers>
#include <optional>
#include <string>
#include <system_error>
#include <tuple>
#include <utility>
#include <vector>

#include "VertexingDataHelper.hpp"

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {

using Covariance = BoundMatrix;
using Propagator = Acts::Propagator<EigenStepper<>>;
using Linearizer = HelicalTrackLinearizer;

// Create a test context
GeometryContext geoContext = GeometryContext::dangerouslyDefaultConstruct();
MagneticFieldContext magFieldContext = MagneticFieldContext();

const std::string toolString = "AMVF";

class FixedVertexSeedFinder final : public IVertexFinder {
 public:
  struct State {};

  explicit FixedVertexSeedFinder(std::vector<Vertex> seeds)
      : m_seeds(std::move(seeds)) {}

  Result<std::vector<Vertex>> find(
      const std::vector<InputTrack>& /*trackVector*/,
      const VertexingOptions& vertexingOptions,
      IVertexFinder::State& /*state*/) const override {
    ++calls;
    lastReference = vertexingOptions.seedReference;
    receivedFitConstraint = vertexingOptions.fitConstraint.has_value();
    return m_seeds;
  }

  IVertexFinder::State makeState(
      const MagneticFieldContext& /*mctx*/) const override {
    return IVertexFinder::State{State{}};
  }

  void setTracksToRemove(
      IVertexFinder::State& /*state*/,
      const std::vector<InputTrack>& /*removedTracks*/) const override {}

  mutable std::size_t calls = 0;
  mutable std::optional<Vertex> lastReference;
  mutable bool receivedFitConstraint = false;

 private:
  std::vector<Vertex> m_seeds;
};

BOOST_AUTO_TEST_SUITE(VertexingSuite)

BOOST_AUTO_TEST_CASE(adaptive_multi_vertex_finder_prior_modes) {
  auto bField = std::make_shared<ConstantBField>(Vector3::Zero());
  auto propagator = std::make_shared<Propagator>(EigenStepper<>(bField));
  ImpactPointEstimator ipEstimator(
      ImpactPointEstimator::Config(bField, propagator));
  Linearizer::Config linearizerConfig;
  linearizerConfig.propagator = propagator;
  Linearizer linearizer(linearizerConfig);

  Vertex seed(Vector3(0.02_mm, 0.03_mm, 0.04_mm));
  seed.setCovariance(SquareMatrix3::Identity() * 1e-4_mm2);
  auto seedFinder =
      std::make_shared<FixedVertexSeedFinder>(std::vector<Vertex>{seed});

  auto makeFinder = [&](bool useSeedConstraint) {
    AdaptiveMultiVertexFitter::Config fitterConfig(ipEstimator);
    fitterConfig.extractParameters.connect<&InputTrack::extractParameters>();
    fitterConfig.trackLinearizer.connect<&Linearizer::linearizeTrack>(
        &linearizer);
    AdaptiveMultiVertexFinder::Config finderConfig(
        AdaptiveMultiVertexFitter(fitterConfig), seedFinder, ipEstimator,
        bField);
    finderConfig.extractParameters.connect<&InputTrack::extractParameters>();
    finderConfig.useSeedConstraint = useSeedConstraint;
    finderConfig.addSingleTrackVertices = true;
    finderConfig.maxIterations = 1;
    return AdaptiveMultiVertexFinder(std::move(finderConfig));
  };

  auto surface = Surface::makeShared<PerigeeSurface>(Vector3::Zero());
  BoundMatrix covariance = BoundMatrix::Zero();
  covariance.diagonal() << 0.01_mm2, 0.01_mm2, 1e-4, 1e-4, 1e-6, 1.;
  std::vector<BoundTrackParameters> tracks;
  for (double phi : {0., std::numbers::pi / 2.}) {
    BoundVector parameters;
    parameters << 0., 0., phi, std::numbers::pi / 2., 1. / 1_GeV, 0.;
    tracks.emplace_back(surface, parameters, covariance,
                        ParticleHypothesis::pion());
  }
  std::vector<InputTrack> inputTracks;
  for (const auto& track : tracks) {
    inputTracks.emplace_back(&track);
  }

  // An arbitrary reference is forwarded to the seeder, while its broad
  // covariance remains independent of the fit's numerical initialization.
  Vertex reference(Vector3(3_mm, -4_mm, 5_mm));
  reference.setCovariance(SquareMatrix3::Identity() * 1e8_mm2);
  VertexingOptions options(geoContext, magFieldContext);
  options.seedReference = reference;
  auto finder = makeFinder(false);
  auto state = finder.makeState(magFieldContext);
  auto freeResult = finder.find(inputTracks, options, state);
  BOOST_REQUIRE(freeResult.ok());
  BOOST_REQUIRE_EQUAL(freeResult->size(), 1u);
  BOOST_REQUIRE(seedFinder->lastReference.has_value());
  CHECK_CLOSE_ABS(seedFinder->lastReference->fullPosition(),
                  reference.fullPosition(), 1e-12);
  CHECK_CLOSE_ABS(seedFinder->lastReference->fullCovariance(),
                  reference.fullCovariance(), 1e-12);
  BOOST_CHECK(!seedFinder->receivedFitConstraint);
  CHECK_CLOSE_ABS(freeResult->front().position(), Vector3::Zero(), 1e-5_mm);

  // Removing only the reference does not alter the fit of the same seed.
  options.seedReference.reset();
  state = finder.makeState(magFieldContext);
  auto noReferenceResult = finder.find(inputTracks, options, state);
  BOOST_REQUIRE(noReferenceResult.ok());
  BOOST_REQUIRE_EQUAL(noReferenceResult->size(), 1u);
  CHECK_CLOSE_ABS(noReferenceResult->front().position(),
                  freeResult->front().position(), 1e-12);
  CHECK_CLOSE_ABS(noReferenceResult->front().covariance(),
                  freeResult->front().covariance(), 1e-12);

  // The external prior is distinct from both reference and seed. It must
  // constrain the reconstructed position without being sent to the seeder.
  Vertex prior(Vector3(0.08_mm, -0.05_mm, 0.1_mm));
  prior.setCovariance(SquareMatrix3::Identity() * 1e-5_mm2);
  options.seedReference = reference;
  options.fitConstraint = prior;
  state = finder.makeState(magFieldContext);
  auto constrainedResult = finder.find(inputTracks, options, state);
  BOOST_REQUIRE(constrainedResult.ok());
  BOOST_REQUIRE_EQUAL(constrainedResult->size(), 1u);
  BOOST_CHECK(!seedFinder->receivedFitConstraint);
  CHECK_CLOSE_ABS(seedFinder->lastReference->fullPosition(),
                  reference.fullPosition(), 1e-12);
  CHECK_CLOSE_ABS(constrainedResult->front().position(), prior.position(),
                  1e-3_mm);
  BOOST_CHECK_LT(constrainedResult->front().covariance().trace(),
                 freeResult->front().covariance().trace());

  options.fitConstraint->setCovariance(SquareMatrix3::Zero());
  auto callsBefore = seedFinder->calls;
  auto invalidPrior = finder.find(inputTracks, options, state);
  BOOST_REQUIRE(!invalidPrior.ok());
  BOOST_CHECK(invalidPrior.error() == VertexingError::InvalidInput);
  BOOST_CHECK_EQUAL(seedFinder->calls, callsBefore);
  options.fitConstraint = prior;

  auto seedConstrainedFinder = makeFinder(true);
  state = seedConstrainedFinder.makeState(magFieldContext);
  callsBefore = seedFinder->calls;
  auto conflict = seedConstrainedFinder.find(inputTracks, options, state);
  BOOST_REQUIRE(!conflict.ok());
  BOOST_CHECK(conflict.error() == VertexingError::InvalidInput);
  BOOST_CHECK_EQUAL(seedFinder->calls, callsBefore);

  options.fitConstraint.reset();
  auto seedConstrainedResult =
      seedConstrainedFinder.find(inputTracks, options, state);
  BOOST_REQUIRE(seedConstrainedResult.ok());
  BOOST_REQUIRE_EQUAL(seedConstrainedResult->size(), 1u);
  CHECK_CLOSE_ABS(seedConstrainedResult->front().position(), seed.position(),
                  1e-3_mm);
  BOOST_CHECK_GT(seedConstrainedResult->front().fitQuality().second,
                 freeResult->front().fitQuality().second);

  // Single-track acceptance still requires a physical or seed-derived prior;
  // the broad unconstrained initialization must not enable it.
  inputTracks.pop_back();
  state = finder.makeState(magFieldContext);
  auto freeSingleTrack = finder.find(inputTracks, options, state);
  BOOST_REQUIRE(freeSingleTrack.ok());
  BOOST_CHECK(freeSingleTrack->empty());
  state = seedConstrainedFinder.makeState(magFieldContext);
  auto seedSingleTrack =
      seedConstrainedFinder.find(inputTracks, options, state);
  BOOST_REQUIRE(seedSingleTrack.ok());
  BOOST_CHECK_EQUAL(seedSingleTrack->size(), 1u);
}

/// @brief AMVF test with Gaussian seed finder
BOOST_AUTO_TEST_CASE(adaptive_multi_vertex_finder_test) {
  // Set debug mode
  bool debugMode = false;
  // Set up constant B-Field
  auto bField = std::make_shared<ConstantBField>(Vector3(0., 0., 2_T));

  // Set up EigenStepper
  // EigenStepper<> stepper(bField);
  EigenStepper<> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);

  // IP 3D Estimator
  ImpactPointEstimator::Config ipEstimatorCfg(bField, propagator);
  ImpactPointEstimator ipEstimator(ipEstimatorCfg);

  std::vector<double> temperatures{
      8., 4., 2., std::numbers::sqrt2, std::sqrt(3. / 2.), 1.};
  AnnealingUtility::Config annealingConfig;
  annealingConfig.setOfTemperatures = temperatures;
  AnnealingUtility annealingUtility(annealingConfig);

  using Fitter = AdaptiveMultiVertexFitter;

  Fitter::Config fitterCfg(ipEstimator);

  fitterCfg.annealingTool = annealingUtility;

  // Linearizer for BoundTrackParameters type test
  Linearizer::Config ltConfig;
  ltConfig.propagator = propagator;
  Linearizer linearizer(ltConfig);

  // Test smoothing
  fitterCfg.doSmoothing = true;
  fitterCfg.extractParameters.connect<&InputTrack::extractParameters>();
  fitterCfg.trackLinearizer.connect<&Linearizer::linearizeTrack>(&linearizer);

  Fitter fitter(fitterCfg);

  GaussianTrackDensity::Config densityCfg;
  densityCfg.extractParameters.connect<&InputTrack::extractParameters>();
  auto seedFinder = std::make_shared<TrackDensityVertexFinder>(
      TrackDensityVertexFinder::Config{Acts::GaussianTrackDensity(densityCfg)});

  AdaptiveMultiVertexFinder::Config finderConfig(std::move(fitter), seedFinder,
                                                 ipEstimator, bField);
  finderConfig.extractParameters.connect<&InputTrack::extractParameters>();
  finderConfig.useSeedConstraint = true;

  AdaptiveMultiVertexFinder finder(std::move(finderConfig));
  IVertexFinder::State state = finder.makeState(magFieldContext);

  auto csvData = readTracksAndVertexCSV(toolString);
  std::vector<BoundTrackParameters> tracks = std::get<TracksData>(csvData);

  if (debugMode) {
    std::cout << "Number of tracks in event: " << tracks.size() << std::endl;
    int maxCout = 10;
    int count = 0;
    for (const auto& trk : tracks) {
      std::cout << count << ". track: " << std::endl;
      std::cout << "params: " << trk << std::endl;
      count++;
      if (count == maxCout) {
        break;
      }
    }
  }

  std::vector<InputTrack> inputTracks;
  for (const auto& trk : tracks) {
    inputTracks.emplace_back(&trk);
  }

  // TODO: test without using beam spot constraint
  Vertex bsConstr = std::get<BeamSpotData>(csvData);
  VertexingOptions vertexingOptions(geoContext, magFieldContext);
  vertexingOptions.seedReference = bsConstr;

  auto t1 = std::chrono::system_clock::now();
  auto findResult = finder.find(inputTracks, vertexingOptions, state);
  auto t2 = std::chrono::system_clock::now();

  auto timediff =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  if (!findResult.ok()) {
    std::cout << findResult.error().message() << std::endl;
  }

  BOOST_CHECK(findResult.ok());

  std::vector<Vertex> allVertices = *findResult;

  if (debugMode) {
    std::cout << "Time needed: " << timediff << " ms." << std::endl;
    std::cout << "Number of vertices reconstructed: " << allVertices.size()
              << std::endl;

    int count = 0;
    for (const auto& vtx : allVertices) {
      count++;
      std::cout << count << ". Vertex at position: " << vtx.position()[0]
                << ", " << vtx.position()[1] << ", " << vtx.position()[2]
                << std::endl;
      std::cout << count << ". Vertex with cov: " << vtx.covariance()
                << std::endl;
      std::cout << "\t with n tracks: " << vtx.tracks().size() << std::endl;
    }
  }

  // Test expected outcomes from athena implementation
  // Number of reconstructed vertices
  auto verticesInfo = std::get<VerticesData>(csvData);
  const int expNRecoVertices = verticesInfo.size();

  BOOST_CHECK_EQUAL(allVertices.size(), expNRecoVertices);

  double relTol = 1e-2;
  double small = 1e-3;
  for (int i = 0; i < expNRecoVertices; i++) {
    auto recoVtx = allVertices[i];
    auto expVtx = verticesInfo[i];
    CHECK_CLOSE_OR_SMALL(recoVtx.position(), expVtx.position, relTol, small);
    CHECK_CLOSE_OR_SMALL(recoVtx.covariance(), expVtx.covariance, relTol,
                         small);
    BOOST_CHECK_EQUAL(recoVtx.tracks().size(), expVtx.nTracks);
    CHECK_CLOSE_OR_SMALL(recoVtx.tracks()[0].trackWeight, expVtx.trk1Weight,
                         relTol, small);
    CHECK_CLOSE_OR_SMALL(recoVtx.tracks()[0].vertexCompatibility,
                         expVtx.trk1Comp, relTol, small);
  }
}

// Dummy user-defined InputTrackStub type
struct InputTrackStub {
  InputTrackStub(const BoundTrackParameters& params, int id)
      : m_parameters(params), m_id(id) {}

  const BoundTrackParameters& parameters() const { return m_parameters; }
  // store e.g. link to original objects here

  int id() const { return m_id; }

 private:
  BoundTrackParameters m_parameters;

  // Some test track ID
  int m_id;
};

/// @brief AMVF test with user-defined input track type
BOOST_AUTO_TEST_CASE(adaptive_multi_vertex_finder_usertype_test) {
  // Set debug mode
  bool debugMode = false;
  // Set up constant B-Field
  auto bField = std::make_shared<ConstantBField>(Vector3(0., 0., 2_T));

  // Set up EigenStepper
  // EigenStepper<> stepper(bField);
  EigenStepper<> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);

  // Create a custom std::function to extract BoundTrackParameters from
  // user-defined InputTrackStub
  auto extractParameters = [](const InputTrack& track) {
    return track.as<InputTrackStub>()->parameters();
  };

  // IP 3D Estimator
  ImpactPointEstimator::Config ipEstimatorCfg(bField, propagator);
  ImpactPointEstimator ipEstimator(ipEstimatorCfg);

  std::vector<double> temperatures{
      8., 4., 2., std::numbers::sqrt2, std::sqrt(3. / 2.), 1.};
  AnnealingUtility::Config annealingConfig;
  annealingConfig.setOfTemperatures = temperatures;
  AnnealingUtility annealingUtility(annealingConfig);

  using Fitter = AdaptiveMultiVertexFitter;

  Fitter::Config fitterCfg(ipEstimator);

  fitterCfg.annealingTool = annealingUtility;

  // Linearizer
  Linearizer::Config ltConfig;
  ltConfig.propagator = propagator;
  Linearizer linearizer(ltConfig);

  // Test smoothing
  fitterCfg.doSmoothing = true;
  fitterCfg.extractParameters.connect(extractParameters);
  fitterCfg.trackLinearizer.connect<&Linearizer::linearizeTrack>(&linearizer);

  Fitter fitter(fitterCfg);

  GaussianTrackDensity::Config densityCfg;
  densityCfg.extractParameters.connect(extractParameters);
  auto seedFinder = std::make_shared<TrackDensityVertexFinder>(
      TrackDensityVertexFinder::Config{Acts::GaussianTrackDensity(densityCfg)});

  AdaptiveMultiVertexFinder::Config finderConfig(
      std::move(fitter), std::move(seedFinder), ipEstimator, bField);
  finderConfig.extractParameters.connect(extractParameters);
  finderConfig.useSeedConstraint = true;

  AdaptiveMultiVertexFinder finder(std::move(finderConfig));
  IVertexFinder::State state = finder.makeState(magFieldContext);

  auto csvData = readTracksAndVertexCSV(toolString);
  auto tracks = std::get<TracksData>(csvData);

  std::vector<InputTrackStub> userTracks;
  int idCount = 0;
  for (const auto& trk : tracks) {
    userTracks.push_back(InputTrackStub(trk, idCount));
    idCount++;
  }

  if (debugMode) {
    std::cout << "Number of tracks in event: " << tracks.size() << std::endl;
    int maxCout = 10;
    int count = 0;
    for (const auto& trk : tracks) {
      std::cout << count << ". track: " << std::endl;
      std::cout << "params: " << trk << std::endl;
      count++;
      if (count == maxCout) {
        break;
      }
    }
  }

  std::vector<InputTrack> userInputTracks;
  for (const auto& trk : userTracks) {
    userInputTracks.emplace_back(&trk);
  }

  Vertex constraintVtx;
  constraintVtx.setPosition(std::get<BeamSpotData>(csvData).position());
  constraintVtx.setCovariance(std::get<BeamSpotData>(csvData).covariance());

  VertexingOptions vertexingOptions(geoContext, magFieldContext);
  vertexingOptions.seedReference = constraintVtx;

  auto findResult = finder.find(userInputTracks, vertexingOptions, state);

  if (!findResult.ok()) {
    std::cout << findResult.error().message() << std::endl;
  }

  BOOST_CHECK(findResult.ok());

  std::vector<Vertex> allVertices = *findResult;

  if (debugMode) {
    std::cout << "Number of vertices reconstructed: " << allVertices.size()
              << std::endl;

    int count = 0;
    for (const auto& vtx : allVertices) {
      count++;
      std::cout << count << ". Vertex at position: " << vtx.position()[0]
                << ", " << vtx.position()[1] << ", " << vtx.position()[2]
                << std::endl;
      std::cout << count << ". Vertex with cov: " << vtx.covariance()
                << std::endl;
      std::cout << "\t with n tracks: " << vtx.tracks().size() << std::endl;
    }
    for (auto& trk : allVertices[0].tracks()) {
      std::cout << "Track ID at first vertex: "
                << trk.originalParams.as<InputTrackStub>()->id() << std::endl;
    }
  }

  auto verticesInfo = std::get<VerticesData>(csvData);
  const int expNRecoVertices = verticesInfo.size();

  BOOST_CHECK_EQUAL(allVertices.size(), expNRecoVertices);

  double relTol = 1e-2;
  double small = 1e-3;
  for (int i = 0; i < expNRecoVertices; i++) {
    auto recoVtx = allVertices[i];
    auto expVtx = verticesInfo[i];
    CHECK_CLOSE_OR_SMALL(recoVtx.position(), expVtx.position, relTol, small);
    CHECK_CLOSE_OR_SMALL(recoVtx.covariance(), expVtx.covariance, relTol,
                         small);
    BOOST_CHECK_EQUAL(recoVtx.tracks().size(), expVtx.nTracks);
    CHECK_CLOSE_OR_SMALL(recoVtx.tracks()[0].trackWeight, expVtx.trk1Weight,
                         relTol, small);
    CHECK_CLOSE_OR_SMALL(recoVtx.tracks()[0].vertexCompatibility,
                         expVtx.trk1Comp, relTol, small);
  }
}

/// @brief AMVF test with grid seed finder
BOOST_AUTO_TEST_CASE(adaptive_multi_vertex_finder_grid_seed_finder_test) {
  // Set debug mode
  bool debugMode = false;
  if (debugMode) {
    std::cout << "Starting AMVF test with grid seed finder..." << std::endl;
  }
  // Set up constant B-Field
  auto bField = std::make_shared<ConstantBField>(Vector3(0., 0., 2_T));

  // Set up EigenStepper
  // EigenStepper<> stepper(bField);
  EigenStepper<> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);

  // IP Estimator
  ImpactPointEstimator::Config ipEstCfg(bField, propagator);
  ImpactPointEstimator ipEst(ipEstCfg);

  std::vector<double> temperatures{
      8., 4., 2., std::numbers::sqrt2, std::sqrt(3. / 2.), 1.};
  AnnealingUtility::Config annealingConfig;
  annealingConfig.setOfTemperatures = temperatures;
  AnnealingUtility annealingUtility(annealingConfig);

  using Fitter = AdaptiveMultiVertexFitter;

  Fitter::Config fitterCfg(ipEst);

  fitterCfg.annealingTool = annealingUtility;

  // Linearizer for BoundTrackParameters type test
  Linearizer::Config ltConfig;
  ltConfig.propagator = propagator;
  Linearizer linearizer(ltConfig);

  // Test smoothing
  fitterCfg.doSmoothing = true;
  fitterCfg.extractParameters.connect<&InputTrack::extractParameters>();
  fitterCfg.trackLinearizer.connect<&Linearizer::linearizeTrack>(&linearizer);

  Fitter fitter(fitterCfg);

  using SeedFinder = GridDensityVertexFinder;
  GaussianGridTrackDensity::Config gDensityConfig(250, 4000, 55);
  GaussianGridTrackDensity gDensity(gDensityConfig);
  SeedFinder::Config seedFinderCfg(gDensity);
  seedFinderCfg.cacheGridStateForTrackRemoval = true;
  seedFinderCfg.extractParameters.connect<&InputTrack::extractParameters>();

  auto seedFinder = std::make_shared<SeedFinder>(seedFinderCfg);

  AdaptiveMultiVertexFinder::Config finderConfig(
      std::move(fitter), std::move(seedFinder), ipEst, bField);
  finderConfig.extractParameters.connect<&InputTrack::extractParameters>();
  finderConfig.useSeedConstraint = true;

  AdaptiveMultiVertexFinder finder(std::move(finderConfig));
  IVertexFinder::State state = finder.makeState(magFieldContext);

  auto csvData = readTracksAndVertexCSV(toolString);
  auto tracks = std::get<TracksData>(csvData);

  if (debugMode) {
    std::cout << "Number of tracks in event: " << tracks.size() << std::endl;
    int maxCout = 10;
    int count = 0;
    for (const auto& trk : tracks) {
      std::cout << count << ". track: " << std::endl;
      std::cout << "params: " << trk << std::endl;
      count++;
      if (count == maxCout) {
        break;
      }
    }
  }

  std::vector<InputTrack> inputTracks;
  for (const auto& trk : tracks) {
    inputTracks.emplace_back(&trk);
  }

  // TODO: test using beam spot constraint
  Vertex bsConstr = std::get<BeamSpotData>(csvData);
  VertexingOptions vertexingOptions(geoContext, magFieldContext);
  vertexingOptions.seedReference = bsConstr;

  auto t1 = std::chrono::system_clock::now();
  auto findResult = finder.find(inputTracks, vertexingOptions, state);
  auto t2 = std::chrono::system_clock::now();

  auto timediff =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  if (!findResult.ok()) {
    std::cout << findResult.error().message() << std::endl;
  }

  BOOST_CHECK(findResult.ok());

  std::vector<Vertex> allVertices = *findResult;

  if (debugMode) {
    std::cout << "Time needed: " << timediff << " ms." << std::endl;
    std::cout << "Number of vertices reconstructed: " << allVertices.size()
              << std::endl;

    int count = 0;
    for (const auto& vtx : allVertices) {
      count++;
      std::cout << count << ". Vertex at position: " << vtx.position()[0]
                << ", " << vtx.position()[1] << ", " << vtx.position()[2]
                << std::endl;
      std::cout << count << ". Vertex with cov: " << vtx.covariance()
                << std::endl;
      std::cout << "\t with n tracks: " << vtx.tracks().size() << std::endl;
    }
  }
  // Test expected outcomes from athena implementation
  // Number of reconstructed vertices
  auto verticesInfo = std::get<VerticesData>(csvData);
  const int expNRecoVertices = verticesInfo.size();

  BOOST_CHECK_EQUAL(allVertices.size(), expNRecoVertices);
  std::vector<bool> vtxFound(expNRecoVertices, false);

  for (const auto& vtx : allVertices) {
    double vtxZ = vtx.position()[2];
    double diffZ = 1e5;
    int foundVtxIdx = -1;
    for (int i = 0; i < expNRecoVertices; i++) {
      if (!vtxFound[i]) {
        if (std::abs(vtxZ - verticesInfo[i].position[2]) < diffZ) {
          diffZ = std::abs(vtxZ - verticesInfo[i].position[2]);
          foundVtxIdx = i;
        }
      }
    }
    if (diffZ < 0.5_mm) {
      vtxFound[foundVtxIdx] = true;
      CHECK_CLOSE_ABS(vtx.tracks().size(), verticesInfo[foundVtxIdx].nTracks,
                      1);
    }
  }
  for (bool found : vtxFound) {
    BOOST_CHECK_EQUAL(found, true);
  }
}

/// @brief AMVF test with adaptive grid seed finder
BOOST_AUTO_TEST_CASE(
    adaptive_multi_vertex_finder_adaptive_grid_seed_finder_test) {
  // Set debug mode
  bool debugMode = false;
  if (debugMode) {
    std::cout << "Starting AMVF test with adaptive grid seed finder..."
              << std::endl;
  }
  // Set up constant B-Field
  auto bField = std::make_shared<ConstantBField>(Vector3(0., 0., 2_T));

  // Set up EigenStepper
  // EigenStepper<> stepper(bField);
  EigenStepper<> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);

  // IP Estimator
  ImpactPointEstimator::Config ipEstCfg(bField, propagator);
  ImpactPointEstimator ipEst(ipEstCfg);

  std::vector<double> temperatures{
      8., 4., 2., std::numbers::sqrt2, std::sqrt(3. / 2.), 1.};
  AnnealingUtility::Config annealingConfig;
  annealingConfig.setOfTemperatures = temperatures;
  AnnealingUtility annealingUtility(annealingConfig);

  using Fitter = AdaptiveMultiVertexFitter;

  Fitter::Config fitterCfg(ipEst);

  fitterCfg.annealingTool = annealingUtility;

  // Linearizer for BoundTrackParameters type test
  Linearizer::Config ltConfig;
  ltConfig.propagator = propagator;
  Linearizer linearizer(ltConfig);

  // Test smoothing
  fitterCfg.doSmoothing = true;
  fitterCfg.extractParameters.connect<&InputTrack::extractParameters>();
  fitterCfg.trackLinearizer.connect<&Linearizer::linearizeTrack>(&linearizer);

  Fitter fitter(fitterCfg);

  // Grid density used during vertex seed finding
  AdaptiveGridTrackDensity::Config gridDensityCfg;
  // force track to have exactly spatialTrkGridSize spatial bins for testing
  // purposes
  gridDensityCfg.spatialTrkGridSizeRange = {55, 55};
  gridDensityCfg.spatialBinExtent = 0.05;
  AdaptiveGridTrackDensity gridDensity(gridDensityCfg);

  using SeedFinder = AdaptiveGridDensityVertexFinder;
  SeedFinder::Config seedFinderCfg(gridDensity);
  seedFinderCfg.cacheGridStateForTrackRemoval = true;
  seedFinderCfg.extractParameters.connect<&InputTrack::extractParameters>();

  auto seedFinder = std::make_shared<SeedFinder>(seedFinderCfg);

  AdaptiveMultiVertexFinder::Config finderConfig(
      std::move(fitter), std::move(seedFinder), ipEst, bField);
  finderConfig.extractParameters.connect<&InputTrack::extractParameters>();
  finderConfig.useSeedConstraint = true;

  AdaptiveMultiVertexFinder finder(std::move(finderConfig));
  IVertexFinder::State state = finder.makeState(magFieldContext);

  auto csvData = readTracksAndVertexCSV(toolString);
  auto tracks = std::get<TracksData>(csvData);

  if (debugMode) {
    std::cout << "Number of tracks in event: " << tracks.size() << std::endl;
    int maxCout = 10;
    int count = 0;
    for (const auto& trk : tracks) {
      std::cout << count << ". track: " << std::endl;
      std::cout << "params: " << trk << std::endl;
      count++;
      if (count == maxCout) {
        break;
      }
    }
  }

  std::vector<InputTrack> inputTracks;
  for (const auto& trk : tracks) {
    inputTracks.emplace_back(&trk);
  }

  Vertex bsConstr = std::get<BeamSpotData>(csvData);
  VertexingOptions vertexingOptions(geoContext, magFieldContext);
  vertexingOptions.seedReference = bsConstr;

  auto t1 = std::chrono::system_clock::now();
  auto findResult = finder.find(inputTracks, vertexingOptions, state);
  auto t2 = std::chrono::system_clock::now();

  auto timediff =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  if (!findResult.ok()) {
    std::cout << findResult.error().message() << std::endl;
  }

  BOOST_CHECK(findResult.ok());

  std::vector<Vertex> allVertices = *findResult;

  if (debugMode) {
    std::cout << "Time needed: " << timediff << " ms." << std::endl;
    std::cout << "Number of vertices reconstructed: " << allVertices.size()
              << std::endl;

    int count = 0;
    for (const auto& vtx : allVertices) {
      count++;
      std::cout << count << ". Vertex at position: " << vtx.position()[0]
                << ", " << vtx.position()[1] << ", " << vtx.position()[2]
                << std::endl;
      std::cout << count << ". Vertex with cov: " << vtx.covariance()
                << std::endl;
      std::cout << "\t with n tracks: " << vtx.tracks().size() << std::endl;
    }
  }
  // Test expected outcomes from athena implementation
  // Number of reconstructed vertices
  auto verticesInfo = std::get<VerticesData>(csvData);
  const int expNRecoVertices = verticesInfo.size();

  BOOST_CHECK_EQUAL(allVertices.size(), expNRecoVertices);
  std::vector<bool> vtxFound(expNRecoVertices, false);

  for (const auto& vtx : allVertices) {
    double vtxZ = vtx.position()[2];
    double diffZ = 1e5;
    int foundVtxIdx = -1;
    for (int i = 0; i < expNRecoVertices; i++) {
      if (!vtxFound[i]) {
        if (std::abs(vtxZ - verticesInfo[i].position[2]) < diffZ) {
          diffZ = std::abs(vtxZ - verticesInfo[i].position[2]);
          foundVtxIdx = i;
        }
      }
    }
    if (diffZ < 0.5_mm) {
      vtxFound[foundVtxIdx] = true;
      CHECK_CLOSE_ABS(vtx.tracks().size(), verticesInfo[foundVtxIdx].nTracks,
                      2);
    }
  }
  for (bool found : vtxFound) {
    BOOST_CHECK_EQUAL(found, true);
  }
}

BOOST_AUTO_TEST_CASE(adaptive_multi_vertex_finder_removes_rejected_seed) {
  auto bField = std::make_shared<ConstantBField>(Vector3(0., 0., 2_T));

  EigenStepper<> stepper(bField);
  auto propagator = std::make_shared<Propagator>(stepper);

  ImpactPointEstimator::Config ipEstCfg(bField, propagator);
  ImpactPointEstimator ipEst(ipEstCfg);

  std::vector<double> temperatures{
      8., 4., 2., std::numbers::sqrt2, std::sqrt(3. / 2.), 1.};
  AnnealingUtility::Config annealingConfig;
  annealingConfig.setOfTemperatures = temperatures;
  AnnealingUtility annealingUtility(annealingConfig);

  using Fitter = AdaptiveMultiVertexFitter;

  Fitter::Config fitterCfg(ipEst);
  fitterCfg.annealingTool = annealingUtility;

  Linearizer::Config ltConfig;
  ltConfig.propagator = propagator;
  Linearizer linearizer(ltConfig);

  fitterCfg.doSmoothing = true;
  fitterCfg.extractParameters.connect<&InputTrack::extractParameters>();
  fitterCfg.trackLinearizer.connect<&Linearizer::linearizeTrack>(&linearizer);

  Fitter fitter(fitterCfg);

  auto csvData = readTracksAndVertexCSV(toolString);
  auto tracks = std::get<TracksData>(csvData);
  auto verticesInfo = std::get<VerticesData>(csvData);
  BOOST_REQUIRE_GE(verticesInfo.size(), 2u);

  std::vector<Vertex> seeds;
  seeds.emplace_back(verticesInfo[0].position);
  seeds.back().setCovariance(verticesInfo[0].covariance);
  seeds.emplace_back(verticesInfo[1].position);
  seeds.back().setCovariance(verticesInfo[1].covariance);

  auto seedFinder = std::make_shared<FixedVertexSeedFinder>(std::move(seeds));

  AdaptiveMultiVertexFinder::Config finderConfig(
      std::move(fitter), std::move(seedFinder), ipEst, bField);
  finderConfig.extractParameters.connect<&InputTrack::extractParameters>();
  finderConfig.useSeedConstraint = true;
  finderConfig.maxMergeVertexSignificance = 1e9;
  finderConfig.maxIterations = 1;

  AdaptiveMultiVertexFinder finder(std::move(finderConfig));
  IVertexFinder::State state = finder.makeState(magFieldContext);

  std::vector<InputTrack> inputTracks;
  for (const auto& trk : tracks) {
    inputTracks.emplace_back(&trk);
  }

  VertexingOptions vertexingOptions(geoContext, magFieldContext);
  vertexingOptions.seedReference = std::get<BeamSpotData>(csvData);

  auto findResult = finder.find(inputTracks, vertexingOptions, state);

  if (!findResult.ok()) {
    std::cout << findResult.error().message() << std::endl;
  }

  BOOST_REQUIRE(findResult.ok());

  std::vector<Vertex> allVertices = *findResult;
  BOOST_CHECK_EQUAL(allVertices.size(), 1u);
  BOOST_CHECK_GT(allVertices.front().tracks().size(), 1u);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

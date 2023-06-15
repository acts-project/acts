// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SingleBoundTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/AnnealingUtility.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/AdaptiveGridDensityVertexFinder.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFinder.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFitter.hpp"
#include "Acts/Vertexing/GaussianTrackDensity.hpp"
#include "Acts/Vertexing/GridDensityVertexFinder.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/TrackDensityVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <system_error>
#include <tuple>
#include <utility>
#include <vector>

#include "VertexingDataHelper.hpp"

namespace Acts {
namespace Test {

using namespace Acts::UnitLiterals;

using Covariance = BoundSymMatrix;
using Propagator = Acts::Propagator<EigenStepper<>>;
using Linearizer = HelicalTrackLinearizer<Propagator>;

// Create a test context
GeometryContext geoContext = GeometryContext();
MagneticFieldContext magFieldContext = MagneticFieldContext();

const std::string toolString = "AMVF";

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
  using IPEstimator = ImpactPointEstimator<BoundTrackParameters, Propagator>;

  IPEstimator::Config ipEstimatorCfg(bField, propagator);
  IPEstimator ipEstimator(ipEstimatorCfg);

  std::vector<double> temperatures{8.0, 4.0, 2.0, 1.4142136, 1.2247449, 1.0};
  AnnealingUtility::Config annealingConfig;
  annealingConfig.setOfTemperatures = temperatures;
  AnnealingUtility annealingUtility(annealingConfig);

  using Fitter = AdaptiveMultiVertexFitter<BoundTrackParameters, Linearizer>;

  Fitter::Config fitterCfg(ipEstimator);

  fitterCfg.annealingTool = annealingUtility;

  // Linearizer for BoundTrackParameters type test
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  // Test smoothing
  fitterCfg.doSmoothing = true;

  Fitter fitter(fitterCfg);

  using SeedFinder =
      TrackDensityVertexFinder<Fitter,
                               GaussianTrackDensity<BoundTrackParameters>>;

  SeedFinder seedFinder;

  using Finder = AdaptiveMultiVertexFinder<Fitter, SeedFinder>;

  Finder::Config finderConfig(std::move(fitter), seedFinder, ipEstimator,
                              std::move(linearizer), bField);

  // TODO: test this as well!
  // finderConfig.useBeamSpotConstraint = false;

  Finder finder(finderConfig);
  Finder::State state;

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

  std::vector<const BoundTrackParameters*> tracksPtr;
  for (const auto& trk : tracks) {
    tracksPtr.push_back(&trk);
  }

  VertexingOptions<BoundTrackParameters> vertexingOptions(geoContext,
                                                          magFieldContext);

  vertexingOptions.vertexConstraint = std::get<BeamSpotData>(csvData);

  auto t1 = std::chrono::system_clock::now();
  auto findResult = finder.find(tracksPtr, vertexingOptions, state);
  auto t2 = std::chrono::system_clock::now();

  auto timediff =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  if (!findResult.ok()) {
    std::cout << findResult.error().message() << std::endl;
  }

  BOOST_CHECK(findResult.ok());

  std::vector<Vertex<BoundTrackParameters>> allVertices = *findResult;

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

  for (int i = 0; i < expNRecoVertices; i++) {
    auto recoVtx = allVertices[i];
    auto expVtx = verticesInfo[i];
    CHECK_CLOSE_ABS(recoVtx.position(), expVtx.position, 0.001_mm);
    CHECK_CLOSE_ABS(recoVtx.covariance(), expVtx.covariance, 0.001_mm);
    BOOST_CHECK_EQUAL(recoVtx.tracks().size(), expVtx.nTracks);
    CHECK_CLOSE_ABS(recoVtx.tracks()[0].trackWeight, expVtx.trk1Weight, 0.003);
    CHECK_CLOSE_ABS(recoVtx.tracks()[0].vertexCompatibility, expVtx.trk1Comp,
                    0.003);
  }
}

// Dummy user-defined InputTrack type
struct InputTrack {
  InputTrack(const BoundTrackParameters& params, int id)
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
  // user-defined InputTrack
  std::function<BoundTrackParameters(InputTrack)> extractParameters =
      [](const InputTrack& params) { return params.parameters(); };

  // IP 3D Estimator
  using IPEstimator = ImpactPointEstimator<InputTrack, Propagator>;

  IPEstimator::Config ipEstimatorCfg(bField, propagator);
  IPEstimator ipEstimator(ipEstimatorCfg);

  std::vector<double> temperatures{8.0, 4.0, 2.0, 1.4142136, 1.2247449, 1.0};
  AnnealingUtility::Config annealingConfig;
  annealingConfig.setOfTemperatures = temperatures;
  AnnealingUtility annealingUtility(annealingConfig);

  using Fitter = AdaptiveMultiVertexFitter<InputTrack, Linearizer>;

  Fitter::Config fitterCfg(ipEstimator);

  fitterCfg.annealingTool = annealingUtility;

  // Linearizer
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  // Test smoothing
  fitterCfg.doSmoothing = true;

  Fitter fitter(fitterCfg, extractParameters);

  using SeedFinder =
      TrackDensityVertexFinder<Fitter, GaussianTrackDensity<InputTrack>>;

  SeedFinder seedFinder(extractParameters);

  using Finder = AdaptiveMultiVertexFinder<Fitter, SeedFinder>;

  Finder::Config finderConfig(std::move(fitter), seedFinder, ipEstimator,
                              std::move(linearizer), bField);
  Finder::State state;

  Finder finder(finderConfig, extractParameters);

  auto csvData = readTracksAndVertexCSV(toolString);
  auto tracks = std::get<TracksData>(csvData);

  std::vector<InputTrack> userTracks;
  int idCount = 0;
  for (const auto& trk : tracks) {
    userTracks.push_back(InputTrack(trk, idCount));
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

  std::vector<const InputTrack*> userTracksPtr;
  for (const auto& trk : userTracks) {
    userTracksPtr.push_back(&trk);
  }

  VertexingOptions<InputTrack> vertexingOptions(geoContext, magFieldContext);

  Vertex<InputTrack> constraintVtx;
  constraintVtx.setPosition(std::get<BeamSpotData>(csvData).position());
  constraintVtx.setCovariance(std::get<BeamSpotData>(csvData).covariance());

  vertexingOptions.vertexConstraint = constraintVtx;

  auto findResult = finder.find(userTracksPtr, vertexingOptions, state);

  if (!findResult.ok()) {
    std::cout << findResult.error().message() << std::endl;
  }

  BOOST_CHECK(findResult.ok());

  std::vector<Vertex<InputTrack>> allVertices = *findResult;

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
      std::cout << "Track ID at first vertex: " << trk.originalParams->id()
                << std::endl;
    }
  }

  auto verticesInfo = std::get<VerticesData>(csvData);
  const int expNRecoVertices = verticesInfo.size();

  BOOST_CHECK_EQUAL(allVertices.size(), expNRecoVertices);

  for (int i = 0; i < expNRecoVertices; i++) {
    auto recoVtx = allVertices[i];
    auto expVtx = verticesInfo[i];
    CHECK_CLOSE_ABS(recoVtx.position(), expVtx.position, 0.001_mm);
    CHECK_CLOSE_ABS(recoVtx.covariance(), expVtx.covariance, 0.001_mm);
    BOOST_CHECK_EQUAL(recoVtx.tracks().size(), expVtx.nTracks);
    CHECK_CLOSE_ABS(recoVtx.tracks()[0].trackWeight, expVtx.trk1Weight, 0.003);
    CHECK_CLOSE_ABS(recoVtx.tracks()[0].vertexCompatibility, expVtx.trk1Comp,
                    0.003);
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
  using IPEstimator = ImpactPointEstimator<BoundTrackParameters, Propagator>;

  IPEstimator::Config ipEstCfg(bField, propagator);
  IPEstimator ipEst(ipEstCfg);

  std::vector<double> temperatures{8.0, 4.0, 2.0, 1.4142136, 1.2247449, 1.0};
  AnnealingUtility::Config annealingConfig;
  annealingConfig.setOfTemperatures = temperatures;
  AnnealingUtility annealingUtility(annealingConfig);

  using Fitter = AdaptiveMultiVertexFitter<BoundTrackParameters, Linearizer>;

  Fitter::Config fitterCfg(ipEst);

  fitterCfg.annealingTool = annealingUtility;

  // Linearizer for BoundTrackParameters type test
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  // Test smoothing
  fitterCfg.doSmoothing = true;

  Fitter fitter(fitterCfg);

  using SeedFinder = GridDensityVertexFinder<4000, 55>;
  SeedFinder::Config seedFinderCfg(250);
  seedFinderCfg.cacheGridStateForTrackRemoval = true;

  SeedFinder seedFinder(seedFinderCfg);

  using Finder = AdaptiveMultiVertexFinder<Fitter, SeedFinder>;

  Finder::Config finderConfig(std::move(fitter), seedFinder, ipEst,
                              std::move(linearizer), bField);

  // TODO: test this as well!
  // finderConfig.useBeamSpotConstraint = false;

  Finder finder(finderConfig);
  Finder::State state;

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

  std::vector<const BoundTrackParameters*> tracksPtr;
  for (const auto& trk : tracks) {
    tracksPtr.push_back(&trk);
  }

  VertexingOptions<BoundTrackParameters> vertexingOptions(geoContext,
                                                          magFieldContext);

  vertexingOptions.vertexConstraint = std::get<BeamSpotData>(csvData);

  auto t1 = std::chrono::system_clock::now();
  auto findResult = finder.find(tracksPtr, vertexingOptions, state);
  auto t2 = std::chrono::system_clock::now();

  auto timediff =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  if (!findResult.ok()) {
    std::cout << findResult.error().message() << std::endl;
  }

  BOOST_CHECK(findResult.ok());

  std::vector<Vertex<BoundTrackParameters>> allVertices = *findResult;

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
      if (not vtxFound[i]) {
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
  using IPEstimator = ImpactPointEstimator<BoundTrackParameters, Propagator>;

  IPEstimator::Config ipEstCfg(bField, propagator);
  IPEstimator ipEst(ipEstCfg);

  std::vector<double> temperatures{8.0, 4.0, 2.0, 1.4142136, 1.2247449, 1.0};
  AnnealingUtility::Config annealingConfig;
  annealingConfig.setOfTemperatures = temperatures;
  AnnealingUtility annealingUtility(annealingConfig);

  using Fitter = AdaptiveMultiVertexFitter<BoundTrackParameters, Linearizer>;

  Fitter::Config fitterCfg(ipEst);

  fitterCfg.annealingTool = annealingUtility;

  // Linearizer for BoundTrackParameters type test
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  // Test smoothing
  fitterCfg.doSmoothing = true;

  Fitter fitter(fitterCfg);

  using SeedFinder = AdaptiveGridDensityVertexFinder<55>;
  SeedFinder::Config seedFinderCfg(0.05);
  seedFinderCfg.cacheGridStateForTrackRemoval = true;

  SeedFinder seedFinder(seedFinderCfg);

  using Finder = AdaptiveMultiVertexFinder<Fitter, SeedFinder>;

  Finder::Config finderConfig(std::move(fitter), seedFinder, ipEst,
                              std::move(linearizer), bField);

  Finder finder(finderConfig);
  Finder::State state;

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

  std::vector<const BoundTrackParameters*> tracksPtr;
  for (const auto& trk : tracks) {
    tracksPtr.push_back(&trk);
  }

  VertexingOptions<BoundTrackParameters> vertexingOptions(geoContext,
                                                          magFieldContext);

  vertexingOptions.vertexConstraint = std::get<BeamSpotData>(csvData);

  auto t1 = std::chrono::system_clock::now();
  auto findResult = finder.find(tracksPtr, vertexingOptions, state);
  auto t2 = std::chrono::system_clock::now();

  auto timediff =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  if (!findResult.ok()) {
    std::cout << findResult.error().message() << std::endl;
  }

  BOOST_CHECK(findResult.ok());

  std::vector<Vertex<BoundTrackParameters>> allVertices = *findResult;

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
      if (not vtxFound[i]) {
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

}  // namespace Test
}  // namespace Acts

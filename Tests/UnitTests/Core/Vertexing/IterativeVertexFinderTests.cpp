// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/FsmwMode1dFinder.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/IterativeVertexFinder.hpp"
#include "Acts/Vertexing/TrackToVertexIPEstimator.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexFinderConcept.hpp"

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

using Covariance = BoundSymMatrix;
using Propagator = Propagator<EigenStepper<ConstantBField>>;
using Linearizer = HelicalTrackLinearizer<Propagator>;

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();

// Vertex x/y position distribution
std::uniform_real_distribution<> vXYDist(-0.1_mm, 0.1_mm);
// Vertex z position distribution
std::uniform_real_distribution<> vZDist(-20_mm, 20_mm);
// Track d0 distribution
std::uniform_real_distribution<> d0Dist(-0.01_mm, 0.01_mm);
// Track z0 distribution
std::uniform_real_distribution<> z0Dist(-0.2_mm, 0.2_mm);
// Track pT distribution
std::uniform_real_distribution<> pTDist(0.4_GeV, 10_GeV);
// Track phi distribution
std::uniform_real_distribution<> phiDist(-M_PI, M_PI);
// Track theta distribution
std::uniform_real_distribution<> thetaDist(1.0, M_PI - 1.0);
// Track charge helper distribution
std::uniform_real_distribution<> qDist(-1, 1);
// Track IP resolution distribution
std::uniform_real_distribution<> resIPDist(0., 100_um);
// Track angular distribution
std::uniform_real_distribution<> resAngDist(0., 0.1);
// Track q/p resolution distribution
std::uniform_real_distribution<> resQoPDist(-0.01, 0.01);
// Number of vertices per test event distribution
std::uniform_int_distribution<> nVertexDist(1, 6);
// Number of tracks per vertex distribution
std::uniform_int_distribution<> nTracksDist(5, 15);

// Dummy user-defined InputTrack type
struct InputTrack {
  InputTrack(const BoundParameters& params) : m_parameters(params) {}

  const BoundParameters& parameters() const { return m_parameters; }

  // store e.g. link to original objects here

 private:
  BoundParameters m_parameters;
};

///
/// @brief Unit test for IterativeVertexFinder for BoundParameters
///
BOOST_AUTO_TEST_CASE(iterative_finder_test) {
  bool debug = false;

  // Set up RNG
  int mySeed = 31415;
  std::mt19937 gen(mySeed);

  // Number of test events
  unsigned int nEvents = 5;  // = nTest

  for (unsigned int iEvent = 0; iEvent < nEvents; ++iEvent) {
    // Set up constant B-Field
    ConstantBField bField(0.0, 0.0, 1_T);

    // Set up Eigenstepper
    EigenStepper<ConstantBField> stepper(bField);

    // Set up propagator with void navigator
    auto propagator = std::make_shared<Propagator>(stepper);

    PropagatorOptions<> pOptions(tgContext, mfContext);

    // Linearizer for BoundParameters type test
    Linearizer::Config ltConfig(bField, propagator, pOptions);
    Linearizer linearizer(ltConfig);

    using BilloirFitter = FullBilloirVertexFitter<BoundParameters, Linearizer>;

    // Set up Billoir Vertex Fitter
    BilloirFitter::Config vertexFitterCfg;

    BilloirFitter bFitter(vertexFitterCfg);

    // Set up all seed finder related things
    TrackToVertexIPEstimator<BoundParameters, Propagator>::Config ipEstCfg(
        propagator, pOptions);
    TrackToVertexIPEstimator<BoundParameters, Propagator> ipEst(ipEstCfg);

    using ZScanSeedFinder = ZScanVertexFinder<BilloirFitter>;

    static_assert(VertexFinderConcept<ZScanSeedFinder>,
                  "Vertex finder does not fulfill vertex finder concept.");

    ZScanSeedFinder::Config sFcfg(std::move(ipEst));

    ZScanSeedFinder sFinder(std::move(sFcfg));

    // Vertex Finder
    using VertexFinder = IterativeVertexFinder<BilloirFitter, ZScanSeedFinder>;

    static_assert(VertexFinderConcept<VertexFinder>,
                  "Vertex finder does not fulfill vertex finder concept.");

    // IP 3D Estimator
    using ImpactPointEstimator =
        ImpactPoint3dEstimator<BoundParameters, Propagator>;

    ImpactPointEstimator::Config ip3dEstCfg(bField, propagator, pOptions);
    ImpactPointEstimator ip3dEst(ip3dEstCfg);

    VertexFinder::Config cfg(std::move(bFitter), std::move(linearizer),
                             std::move(sFinder), std::move(ip3dEst));

    cfg.reassignTracksAfterFirstFit = true;

    VertexFinder finder(cfg);

    // Vector to be filled with all tracks in current event
    std::vector<std::unique_ptr<const BoundParameters>> tracks;

    std::vector<const BoundParameters*> tracksPtr;

    // Vector to be filled with truth vertices for later comparison
    std::vector<Vertex<BoundParameters>> trueVertices;

    // start creating event with nVertices vertices
    unsigned int nVertices = nVertexDist(gen);
    for (unsigned int iVertex = 0; iVertex < nVertices; ++iVertex) {
      // Number of tracks
      unsigned int nTracks = nTracksDist(gen);

      if (debug) {
        std::cout << "Event " << iEvent << ", Vertex " << iVertex << "/"
                  << nVertices << " with " << nTracks << " tracks."
                  << std::endl;
      }
      // Create perigee surface
      std::shared_ptr<PerigeeSurface> perigeeSurface =
          Surface::makeShared<PerigeeSurface>(Vector3D(0., 0., 0.));

      // Create position of vertex and perigee surface
      double x = vXYDist(gen);
      double y = vXYDist(gen);
      double z = vZDist(gen);

      // True vertex
      Vertex<BoundParameters> trueV(Vector3D(x, y, z));
      std::vector<TrackAtVertex<BoundParameters>> tracksAtTrueVtx;

      // Calculate d0 and z0 corresponding to vertex position
      double d0_v = sqrt(x * x + y * y);
      double z0_v = z;

      // Construct random track emerging from vicinity of vertex position
      // Vector to store track objects used for vertex fit
      for (unsigned int iTrack = 0; iTrack < nTracks; iTrack++) {
        // Construct positive or negative charge randomly
        double q = qDist(gen) < 0 ? -1. : 1.;

        // Construct random track parameters
        BoundVector paramVec;
        double z0track = z0_v + z0Dist(gen);
        paramVec << d0_v + d0Dist(gen), z0track, phiDist(gen), thetaDist(gen),
            q / pTDist(gen), 0.;

        // Resolutions
        double res_d0 = resIPDist(gen);
        double res_z0 = resIPDist(gen);
        double res_ph = resAngDist(gen);
        double res_th = resAngDist(gen);
        double res_qp = resQoPDist(gen);

        // Fill vector of track objects with simple covariance matrix
        Covariance covMat;
        covMat << res_d0 * res_d0, 0., 0., 0., 0., 0., 0., res_z0 * res_z0, 0.,
            0., 0., 0., 0., 0., res_ph * res_ph, 0., 0., 0., 0., 0., 0.,
            res_th * res_th, 0., 0., 0., 0., 0., 0., res_qp * res_qp, 0., 0.,
            0., 0., 0., 0., 1.;
        auto params = BoundParameters(tgContext, std::move(covMat), paramVec,
                                      perigeeSurface);

        tracks.push_back(std::make_unique<BoundParameters>(params));

        TrackAtVertex<BoundParameters> trAtVt(0., params, tracks.back().get());
        tracksAtTrueVtx.push_back(trAtVt);
      }

      trueV.setTracksAtVertex(tracksAtTrueVtx);
      trueVertices.push_back(trueV);

    }  // end loop over vertices

    // shuffle list of tracks
    std::shuffle(std::begin(tracks), std::end(tracks), gen);

    for (const auto& trk : tracks) {
      tracksPtr.push_back(trk.get());
    }

    VertexFinderOptions<BoundParameters> vFinderOptions(tgContext, mfContext);

    // find vertices
    auto res = finder.find(tracksPtr, vFinderOptions);

    BOOST_CHECK(res.ok());

    if (!res.ok()) {
      std::cout << res.error().message() << std::endl;
    }

    // Retrieve vertices found by vertex finder
    auto vertexCollection = *res;

    // check if same amount of vertices has been found with tolerance of 2
    CHECK_CLOSE_ABS(vertexCollection.size(), nVertices, 2);

    if (debug) {
      std::cout << "########## RESULT: ########## Event " << iEvent
                << std::endl;
      std::cout << "Number of true vertices: " << nVertices << std::endl;
      std::cout << "Number of reco vertices: " << vertexCollection.size()
                << std::endl;

      int count = 1;
      std::cout << "----- True vertices -----" << std::endl;
      for (const auto& vertex : trueVertices) {
        Vector3D pos = vertex.position();
        std::cout << count << ". True Vertex:\t Position:"
                  << "(" << pos[eX] << "," << pos[eY] << "," << pos[eZ] << ")"
                  << std::endl;
        std::cout << "Number of tracks: " << vertex.tracks().size() << std::endl
                  << std::endl;
        count++;
      }
      std::cout << "----- Reco vertices -----" << std::endl;
      count = 1;
      for (const auto& vertex : vertexCollection) {
        Vector3D pos = vertex.position();
        std::cout << count << ". Reco Vertex:\t Position:"
                  << "(" << pos[eX] << "," << pos[eY] << "," << pos[eZ] << ")"
                  << std::endl;
        std::cout << "Number of tracks: " << vertex.tracks().size() << std::endl
                  << std::endl;
        count++;
      }
    }

    // Check if all vertices have been found with close z-values
    bool allVerticesFound = true;
    for (const auto& trueVertex : trueVertices) {
      SpacePointVector truePos = trueVertex.fullPosition();
      bool currentVertexFound = false;
      for (const auto& recoVertex : vertexCollection) {
        SpacePointVector recoPos = recoVertex.fullPosition();
        // check only for close z distance
        double zDistance = std::abs(truePos[eZ] - recoPos[eZ]);
        if (zDistance < 2_mm) {
          currentVertexFound = true;
        }
      }
      if (!currentVertexFound) {
        allVerticesFound = false;
      }
    }

    // check if found vertices have compatible z values
    BOOST_TEST(allVerticesFound);
  }
}

///
/// @brief Unit test for IterativeVertexFinder
///        for user defined InputTrack track type
///
BOOST_AUTO_TEST_CASE(iterative_finder_test_user_track_type) {
  bool debug = false;

  // Set up RNG
  int mySeed = 31415;
  std::mt19937 gen(mySeed);

  // Number of test events
  unsigned int nEvents = 5;  // = nTest

  for (unsigned int iEvent = 0; iEvent < nEvents; ++iEvent) {
    // Set up constant B-Field
    ConstantBField bField(0.0, 0.0, 1_T);

    // Set up Eigenstepper
    EigenStepper<ConstantBField> stepper(bField);

    // Set up propagator with void navigator
    auto propagator = std::make_shared<Propagator>(stepper);

    PropagatorOptions<> pOptions(tgContext, mfContext);

    // Linearizer for user defined InputTrack type test
    Linearizer::Config ltConfigUT(bField, propagator, pOptions);
    Linearizer linearizer(ltConfigUT);

    // Set up vertex fitter for user track type
    using BilloirFitter = FullBilloirVertexFitter<InputTrack, Linearizer>;

    // Create a custom std::function to extract BoundParameters from
    // user-defined InputTrack
    std::function<BoundParameters(InputTrack)> extractParameters =
        [](InputTrack params) { return params.parameters(); };

    // Set up Billoir Vertex Fitter
    BilloirFitter::Config vertexFitterCfg;

    BilloirFitter bFitter(vertexFitterCfg, extractParameters);

    // Set up all seed finder related things
    TrackToVertexIPEstimator<InputTrack, Propagator>::Config ipEstCfg(
        propagator, pOptions);
    TrackToVertexIPEstimator<InputTrack, Propagator> ipEst(ipEstCfg);

    using ZScanSeedFinder = ZScanVertexFinder<BilloirFitter>;
    ZScanSeedFinder::Config sFcfg(std::move(ipEst));

    ZScanSeedFinder sFinder(std::move(sFcfg), extractParameters);

    // IP 3D Estimator
    using ImpactPointEstimator = ImpactPoint3dEstimator<InputTrack, Propagator>;

    ImpactPointEstimator::Config ip3dEstCfg(bField, propagator, pOptions);
    ImpactPointEstimator ip3dEst(ip3dEstCfg);

    // Vertex Finder
    using VertexFinder = IterativeVertexFinder<BilloirFitter, ZScanSeedFinder>;
    VertexFinder::Config cfg(std::move(bFitter), std::move(linearizer),
                             std::move(sFinder), std::move(ip3dEst));
    cfg.reassignTracksAfterFirstFit = true;

    VertexFinder finder(cfg, extractParameters);

    // Same for user track type tracks
    std::vector<std::unique_ptr<const InputTrack>> tracks;

    std::vector<const InputTrack*> tracksPtr;

    // Vector to be filled with truth vertices for later comparison
    std::vector<Vertex<InputTrack>> trueVertices;

    // start creating event with nVertices vertices
    unsigned int nVertices = nVertexDist(gen);
    for (unsigned int iVertex = 0; iVertex < nVertices; ++iVertex) {
      // Number of tracks
      unsigned int nTracks = nTracksDist(gen);

      if (debug) {
        std::cout << "Event " << iEvent << ", Vertex " << iVertex << "/"
                  << nVertices << " with " << nTracks << " tracks."
                  << std::endl;
      }
      // Create perigee surface
      std::shared_ptr<PerigeeSurface> perigeeSurface =
          Surface::makeShared<PerigeeSurface>(Vector3D(0., 0., 0.));

      // Create position of vertex and perigee surface
      double x = vXYDist(gen);
      double y = vXYDist(gen);
      double z = vZDist(gen);

      // True vertex
      Vertex<InputTrack> trueV(Vector3D(x, y, z));
      std::vector<TrackAtVertex<InputTrack>> tracksAtTrueVtx;

      // Calculate d0 and z0 corresponding to vertex position
      double d0_v = sqrt(x * x + y * y);
      double z0_v = z;

      // Construct random track emerging from vicinity of vertex position
      // Vector to store track objects used for vertex fit
      for (unsigned int iTrack = 0; iTrack < nTracks; iTrack++) {
        // Construct positive or negative charge randomly
        double q = qDist(gen) < 0 ? -1. : 1.;

        // Construct random track parameters
        BoundVector paramVec;
        double z0track = z0_v + z0Dist(gen);
        paramVec << d0_v + d0Dist(gen), z0track, phiDist(gen), thetaDist(gen),
            q / pTDist(gen), 0.;

        // Resolutions
        double res_d0 = resIPDist(gen);
        double res_z0 = resIPDist(gen);
        double res_ph = resAngDist(gen);
        double res_th = resAngDist(gen);
        double res_qp = resQoPDist(gen);

        // Fill vector of track objects now for user track type
        Covariance covMat;

        covMat << res_d0 * res_d0, 0., 0., 0., 0., 0., 0., res_z0 * res_z0, 0.,
            0., 0., 0., 0., 0., res_ph * res_ph, 0., 0., 0., 0., 0., 0.,
            res_th * res_th, 0., 0., 0., 0., 0., 0., res_qp * res_qp, 0., 0.,
            0., 0., 0., 0., 1.;
        auto paramsUT = InputTrack(BoundParameters(tgContext, std::move(covMat),
                                                   paramVec, perigeeSurface));

        tracks.push_back(std::make_unique<InputTrack>(paramsUT));

        auto params = extractParameters(paramsUT);

        TrackAtVertex<InputTrack> trAtVt(0., params, tracks.back().get());
        tracksAtTrueVtx.push_back(trAtVt);
      }

      trueV.setTracksAtVertex(tracksAtTrueVtx);
      trueVertices.push_back(trueV);

    }  // end loop over vertices

    // shuffle list of tracks
    std::shuffle(std::begin(tracks), std::end(tracks), gen);

    for (const auto& trk : tracks) {
      tracksPtr.push_back(trk.get());
    }

    VertexFinderOptions<InputTrack> vFinderOptionsUT(tgContext, mfContext);

    // find vertices
    auto res = finder.find(tracksPtr, vFinderOptionsUT);

    BOOST_CHECK(res.ok());

    if (!res.ok()) {
      std::cout << res.error().message() << std::endl;
    }

    // Retrieve vertices found by vertex finder
    auto vertexCollectionUT = *res;

    // check if same amount of vertices has been found with tolerance of 2
    CHECK_CLOSE_ABS(vertexCollectionUT.size(), nVertices, 2);

    if (debug) {
      std::cout << "########## RESULT: ########## Event " << iEvent
                << std::endl;
      std::cout << "Number of true vertices: " << nVertices << std::endl;
      std::cout << "Number of reco vertices: " << vertexCollectionUT.size()
                << std::endl;

      int count = 1;
      std::cout << "----- True vertices -----" << std::endl;
      for (const auto& vertex : trueVertices) {
        Vector3D pos = vertex.position();
        std::cout << count << ". True Vertex:\t Position:"
                  << "(" << pos[eX] << "," << pos[eY] << "," << pos[eZ] << ")"
                  << std::endl;
        std::cout << "Number of tracks: " << vertex.tracks().size() << std::endl
                  << std::endl;
        count++;
      }
      std::cout << "----- Reco vertices -----" << std::endl;
      count = 1;
      for (const auto& vertex : vertexCollectionUT) {
        Vector3D pos = vertex.position();
        std::cout << count << ". Reco Vertex:\t Position:"
                  << "(" << pos[eX] << "," << pos[eY] << "," << pos[eZ] << ")"
                  << std::endl;
        std::cout << "Number of tracks: " << vertex.tracks().size() << std::endl
                  << std::endl;
        count++;
      }
    }

    // Check if all vertices have been found with close z-values
    bool allVerticesFound = true;
    for (const auto& trueVertex : trueVertices) {
      SpacePointVector truePos = trueVertex.fullPosition();
      bool currentVertexFound = false;
      for (const auto& recoVertex : vertexCollectionUT) {
        SpacePointVector recoPos = recoVertex.fullPosition();
        // check only for close z distance
        double zDistance = std::abs(truePos[eZ] - recoPos[eZ]);
        if (zDistance < 2_mm) {
          currentVertexFound = true;
        }
      }
      if (!currentVertexFound) {
        allVerticesFound = false;
      }
    }

    // check if found vertices have compatible z values
    BOOST_TEST(allVerticesFound);
  }
}

}  // namespace Test
}  // namespace Acts

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

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Vertexing/FsmwMode1dFinder.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/VertexFinderConcept.hpp"
#include "Acts/Vertexing/ZScanVertexFinder.hpp"

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

using Covariance = BoundSymMatrix;
using Propagator = Acts::Propagator<EigenStepper<>>;
using Linearizer_t = HelicalTrackLinearizer<Propagator>;

// Create a test context
GeometryContext geoContext = GeometryContext();
MagneticFieldContext magFieldContext = MagneticFieldContext();

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

///
/// @brief Unit test for ZScanVertexFinder
///
BOOST_AUTO_TEST_CASE(zscan_finder_test) {
  unsigned int nTests = 50;

  for (unsigned int iTest = 0; iTest < nTests; ++iTest) {
    // Number of tracks
    unsigned int nTracks = 30;

    // Set up RNG
    int mySeed = 31415;
    std::mt19937 gen(mySeed);

    // Set up constant B-Field
    auto bField = std::make_shared<ConstantBField>(Vector3{0.0, 0.0, 1_T});

    // Set up Eigenstepper
    EigenStepper<> stepper(bField);

    // Set up propagator with void navigator
    auto propagator = std::make_shared<Propagator>(stepper);

    using BilloirFitter =
        FullBilloirVertexFitter<BoundTrackParameters, Linearizer_t>;

    // Create perigee surface
    std::shared_ptr<PerigeeSurface> perigeeSurface =
        Surface::makeShared<PerigeeSurface>(Vector3(0., 0., 0.));

    // Create position of vertex and perigee surface
    double x = vXYDist(gen);
    double y = vXYDist(gen);
    double z = vZDist(gen);

    // Calculate d0 and z0 corresponding to vertex position
    double d0_v = std::hypot(x, y);
    double z0_v = z;

    // Start constructing nTracks tracks in the following
    std::vector<BoundTrackParameters> tracks;

    // Construct random track emerging from vicinity of vertex position
    // Vector to store track objects used for vertex fit
    for (unsigned int iTrack = 0; iTrack < nTracks; iTrack++) {
      // Construct positive or negative charge randomly
      double q = qDist(gen) < 0 ? -1. : 1.;

      // Construct random track parameters
      BoundVector paramVec = BoundVector::Zero();
      paramVec[eBoundLoc0] = d0_v + d0Dist(gen);
      paramVec[eBoundLoc1] = z0_v + z0Dist(gen);
      paramVec[eBoundPhi] = phiDist(gen);
      paramVec[eBoundTheta] = thetaDist(gen);
      paramVec[eBoundQOverP] = q / pTDist(gen);

      // Resolutions
      double resD0 = resIPDist(gen);
      double resZ0 = resIPDist(gen);
      double resPh = resAngDist(gen);
      double resTh = resAngDist(gen);
      double resQp = resQoPDist(gen);

      // Fill vector of track objects with simple covariance matrix
      Covariance covMat;

      covMat << resD0 * resD0, 0., 0., 0., 0., 0., 0., resZ0 * resZ0, 0., 0.,
          0., 0., 0., 0., resPh * resPh, 0., 0., 0., 0., 0., 0., resTh * resTh,
          0., 0., 0., 0., 0., 0., resQp * resQp, 0., 0., 0., 0., 0., 0., 1.;

      tracks.emplace_back(perigeeSurface, paramVec, std::move(covMat));
    }

    std::vector<const BoundTrackParameters*> tracksPtr;
    for (const auto& trk : tracks) {
      tracksPtr.push_back(&trk);
    }

    using VertexFinder = ZScanVertexFinder<BilloirFitter>;

    static_assert(VertexFinderConcept<VertexFinder>,
                  "Vertex finder does not fulfill vertex finder concept.");

    // Impact point estimator
    using IPEstimator = ImpactPointEstimator<BoundTrackParameters, Propagator>;

    IPEstimator::Config ipEstimatorCfg(bField, propagator);
    IPEstimator ipEstimator(ipEstimatorCfg);

    VertexFinder::Config cfg(ipEstimator);

    VertexFinder finder(cfg);

    VertexingOptions<BoundTrackParameters> vertexingOptions(geoContext,
                                                            magFieldContext);

    VertexFinder::State state;
    auto res = finder.find(tracksPtr, vertexingOptions, state);

    BOOST_CHECK(res.ok());

    if (!res.ok()) {
      std::cout << res.error().message() << std::endl;
    }

    if (res.ok()) {
      BOOST_CHECK(!(*res).empty());
      Vector3 result = (*res).back().position();
      CHECK_CLOSE_ABS(result[eZ], z, 1_mm);
    }
  }
}

// Dummy user-defined InputTrack type
struct InputTrack {
  InputTrack(const BoundTrackParameters& params) : m_parameters(params) {}

  const BoundTrackParameters& parameters() const { return m_parameters; }

  // store e.g. link to original objects here

 private:
  BoundTrackParameters m_parameters;
};

///
/// @brief Unit test for ZScanVertexFinder with user-defined input track type
///
BOOST_AUTO_TEST_CASE(zscan_finder_usertrack_test) {
  unsigned int nTests = 50;

  for (unsigned int iTest = 0; iTest < nTests; ++iTest) {
    // Number of tracks
    unsigned int nTracks = 30;

    // Set up RNG
    int mySeed = 31415;
    std::mt19937 gen(mySeed);

    // Set up constant B-Field
    auto bField = std::make_shared<ConstantBField>(Vector3{0.0, 0.0, 1_T});

    // Set up Eigenstepper
    EigenStepper<> stepper(bField);

    // Set up propagator with void navigator
    auto propagator = std::make_shared<Propagator>(stepper);

    using BilloirFitter = FullBilloirVertexFitter<InputTrack, Linearizer_t>;

    // Create perigee surface
    std::shared_ptr<PerigeeSurface> perigeeSurface =
        Surface::makeShared<PerigeeSurface>(Vector3(0., 0., 0.));

    // Create position of vertex and perigee surface
    double x = vXYDist(gen);
    double y = vXYDist(gen);
    double z = vZDist(gen);

    // Calculate d0 and z0 corresponding to vertex position
    double d0_v = std::hypot(x, y);
    double z0_v = z;

    // Start constructing nTracks tracks in the following
    std::vector<InputTrack> tracks;

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
      double resD0 = resIPDist(gen);
      double resZ0 = resIPDist(gen);
      double resPh = resAngDist(gen);
      double resTh = resAngDist(gen);
      double resQp = resQoPDist(gen);

      // Fill vector of track objects with simple covariance matrix
      Covariance covMat;

      covMat << resD0 * resD0, 0., 0., 0., 0., 0., 0., resZ0 * resZ0, 0., 0.,
          0., 0., 0., 0., resPh * resPh, 0., 0., 0., 0., 0., 0., resTh * resTh,
          0., 0., 0., 0., 0., 0., resQp * resQp, 0., 0., 0., 0., 0., 0., 1.;

      tracks.emplace_back(
          BoundTrackParameters(perigeeSurface, paramVec, std::move(covMat)));
    }

    std::vector<const InputTrack*> tracksPtr;
    for (const auto& trk : tracks) {
      tracksPtr.push_back(&trk);
    }

    using VertexFinder = ZScanVertexFinder<BilloirFitter>;

    static_assert(VertexFinderConcept<VertexFinder>,
                  "Vertex finder does not fulfill vertex finder concept.");

    // Impact point estimator
    using IPEstimator = ImpactPointEstimator<InputTrack, Propagator>;

    IPEstimator::Config ipEstimatorCfg(bField, propagator);
    IPEstimator ipEstimator(ipEstimatorCfg);

    VertexFinder::Config cfg(ipEstimator);

    // Create a custom std::function to extract BoundTrackParameters from
    // user-defined InputTrack
    std::function<BoundTrackParameters(InputTrack)> extractParameters =
        [](const InputTrack& params) { return params.parameters(); };

    VertexFinder finder(cfg, extractParameters);
    VertexFinder::State state;

    VertexingOptions<InputTrack> vertexingOptions(geoContext, magFieldContext);

    auto res = finder.find(tracksPtr, vertexingOptions, state);

    BOOST_CHECK(res.ok());

    if (!res.ok()) {
      std::cout << res.error().message() << std::endl;
    }

    if (res.ok()) {
      BOOST_CHECK(!(*res).empty());
      Vector3 result = (*res).back().position();
      CHECK_CLOSE_ABS(result[eZ], z, 1_mm);
    }
  }
}

}  // namespace Test
}  // namespace Acts

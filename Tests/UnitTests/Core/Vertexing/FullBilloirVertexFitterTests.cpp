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
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/Vertex.hpp"

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

using Covariance = BoundSymMatrix;
using Linearizer = HelicalTrackLinearizer<Propagator<EigenStepper<>>>;

// Create a test context
GeometryContext geoContext = GeometryContext();
MagneticFieldContext magFieldContext = MagneticFieldContext();

/// @brief Unit test for FullBilloirVertexFitter
///
BOOST_AUTO_TEST_CASE(billoir_vertex_fitter_empty_input_test) {
  // Set up constant B-Field
  auto bField = std::make_shared<ConstantBField>(Vector3{0.0, 0.0, 1_T});

  // Set up Eigenstepper
  EigenStepper<> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator<EigenStepper<>>>(stepper);

  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  // Set up Billoir Vertex Fitter
  using VertexFitter =
      FullBilloirVertexFitter<BoundTrackParameters, Linearizer>;
  VertexFitter::Config vertexFitterCfg;
  VertexFitter billoirFitter(vertexFitterCfg);
  VertexFitter::State state(bField->makeCache(magFieldContext));

  // Constraint for vertex fit
  Vertex<BoundTrackParameters> myConstraint;
  // Some abitrary values
  SymMatrix4 myCovMat = SymMatrix4::Zero();
  myCovMat(0, 0) = 30.;
  myCovMat(1, 1) = 30.;
  myCovMat(2, 2) = 30.;
  myCovMat(3, 3) = 30.;
  myConstraint.setFullCovariance(myCovMat);
  myConstraint.setFullPosition(Vector4(0, 0, 0, 0));

  const std::vector<const BoundTrackParameters*> emptyVector;

  VertexingOptions<BoundTrackParameters> vfOptions(geoContext, magFieldContext,
                                                   myConstraint);

  Vertex<BoundTrackParameters> fittedVertex =
      billoirFitter.fit(emptyVector, linearizer, vfOptions, state).value();

  Vector3 origin(0., 0., 0.);
  BOOST_CHECK_EQUAL(fittedVertex.position(), origin);

  SymMatrix4 zeroMat = SymMatrix4::Zero();
  BOOST_CHECK_EQUAL(fittedVertex.fullCovariance(), zeroMat);

  fittedVertex =
      billoirFitter.fit(emptyVector, linearizer, vfOptions, state).value();

  BOOST_CHECK_EQUAL(fittedVertex.position(), origin);
  BOOST_CHECK_EQUAL(fittedVertex.fullCovariance(), zeroMat);
}

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
std::uniform_real_distribution<> resQoPDist(-0.1, 0.1);
// Number of tracks distritbution
std::uniform_int_distribution<> nTracksDist(3, 10);

///
/// @brief Unit test for FullBilloirVertexFitter
/// with default input track type (= BoundTrackParameters)
///
BOOST_AUTO_TEST_CASE(billoir_vertex_fitter_defaulttrack_test) {
  // Set up RNG
  int mySeed = 31415;
  std::mt19937 gen(mySeed);
  // Set up constant B-Field
  auto bField = std::make_shared<ConstantBField>(Vector3{0.0, 0.0, 1_T});

  // Set up Eigenstepper
  EigenStepper<> stepper(bField);
  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator<EigenStepper<>>>(stepper);

  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  // Number of events
  const int nEvents = 10;
  for (int eventIdx = 0; eventIdx < nEvents; ++eventIdx) {
    // Number of tracks
    unsigned int nTracks = nTracksDist(gen);

    // Set up Billoir Vertex Fitter
    using VertexFitter =
        FullBilloirVertexFitter<BoundTrackParameters, Linearizer>;
    VertexFitter::Config vertexFitterCfg;
    VertexFitter billoirFitter(vertexFitterCfg);
    VertexFitter::State state(bField->makeCache(magFieldContext));
    // Constraint for vertex fit
    Vertex<BoundTrackParameters> myConstraint;
    // Some abitrary values
    SymMatrix4 myCovMat = SymMatrix4::Zero();
    myCovMat(0, 0) = 30.;
    myCovMat(1, 1) = 30.;
    myCovMat(2, 2) = 30.;
    myCovMat(3, 3) = 30.;
    myConstraint.setFullCovariance(myCovMat);
    myConstraint.setFullPosition(Vector4(0, 0, 0, 0));
    VertexingOptions<BoundTrackParameters> vfOptions(geoContext,
                                                     magFieldContext);

    VertexingOptions<BoundTrackParameters> vfOptionsConstr(
        geoContext, magFieldContext, myConstraint);
    // Create position of vertex and perigee surface
    double x = vXYDist(gen);
    double y = vXYDist(gen);
    double z = vZDist(gen);

    Vector3 vertexPosition(x, y, z);
    std::shared_ptr<PerigeeSurface> perigeeSurface =
        Surface::makeShared<PerigeeSurface>(Vector3(0., 0., 0.));
    // Calculate d0 and z0 corresponding to vertex position
    double d0V = sqrt(x * x + y * y);
    double z0V = z;

    // Start constructing nTracks tracks in the following
    std::vector<BoundTrackParameters> tracks;

    // Construct random track emerging from vicinity of vertex position
    // Vector to store track objects used for vertex fit
    for (unsigned int iTrack = 0; iTrack < nTracks; iTrack++) {
      // Construct positive or negative charge randomly
      double q = qDist(gen) < 0 ? -1. : 1.;

      // Construct random track parameters
      BoundVector paramVec;
      paramVec << d0V + d0Dist(gen), z0V + z0Dist(gen), phiDist(gen),
          thetaDist(gen), q / pTDist(gen), 0.;

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
      tracks.push_back(
          BoundTrackParameters(perigeeSurface, paramVec, std::move(covMat)));
    }

    std::vector<const BoundTrackParameters*> tracksPtr;
    for (const auto& trk : tracks) {
      tracksPtr.push_back(&trk);
    }

    // Do the actual fit with 4 tracks without constraint
    Vertex<BoundTrackParameters> fittedVertex =
        billoirFitter.fit(tracksPtr, linearizer, vfOptions, state).value();
    if (!fittedVertex.tracks().empty()) {
      CHECK_CLOSE_ABS(fittedVertex.position(), vertexPosition, 1_mm);
    }
    // Do the fit with a constraint
    Vertex<BoundTrackParameters> fittedVertexConstraint =
        billoirFitter.fit(tracksPtr, linearizer, vfOptionsConstr, state)
            .value();
    if (!fittedVertexConstraint.tracks().empty()) {
      CHECK_CLOSE_ABS(fittedVertexConstraint.position(), vertexPosition, 1_mm);
    }

    std::cout << "Fitting nTracks: " << nTracks << std::endl;
    std::cout << "True Vertex: " << x << ", " << y << ", " << z << std::endl;
    std::cout << "Fitted Vertex: " << fittedVertex.position() << std::endl;
    std::cout << "Fitted constraint Vertex: "
              << fittedVertexConstraint.position() << std::endl;
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
/// @brief Unit test for FullBilloirVertexFitter with user-defined InputTrack
/// type
///
BOOST_AUTO_TEST_CASE(billoir_vertex_fitter_usertrack_test) {
  // Set up RNG
  int mySeed = 31415;
  std::mt19937 gen(mySeed);

  // Set up constant B-Field
  auto bField = std::make_shared<ConstantBField>(Vector3{0.0, 0.0, 1_T});

  // Set up Eigenstepper
  EigenStepper<> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator<EigenStepper<>>>(stepper);

  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  const int nEvents = 10;

  for (int eventIdx = 0; eventIdx < nEvents; ++eventIdx) {
    unsigned int nTracks = nTracksDist(gen);

    // Create a custom std::function to extract BoundTrackParameters from
    // user-defined InputTrack
    std::function<BoundTrackParameters(InputTrack)> extractParameters =
        [](InputTrack params) { return params.parameters(); };

    // Set up Billoir Vertex Fitter
    using VertexFitter = FullBilloirVertexFitter<InputTrack, Linearizer>;
    VertexFitter::Config vertexFitterCfg;
    VertexFitter billoirFitter(vertexFitterCfg, extractParameters);
    VertexFitter::State state(bField->makeCache(magFieldContext));

    // Constraint for vertex fit
    Vertex<InputTrack> myConstraint;
    // Some abitrary values
    SymMatrix4 myCovMat = SymMatrix4::Zero();
    myCovMat(0, 0) = 30.;
    myCovMat(1, 1) = 30.;
    myCovMat(2, 2) = 30.;
    myCovMat(3, 3) = 30.;
    myConstraint.setFullCovariance(myCovMat);
    myConstraint.setFullPosition(Vector4(0, 0, 0, 0));

    VertexingOptions<InputTrack> vfOptions(geoContext, magFieldContext);

    VertexingOptions<InputTrack> vfOptionsConstr(geoContext, magFieldContext,
                                                 myConstraint);

    // Create position of vertex and perigee surface
    double x = vXYDist(gen);
    double y = vXYDist(gen);
    double z = vZDist(gen);

    Vector3 vertexPosition(x, y, z);
    std::shared_ptr<PerigeeSurface> perigeeSurface =
        Surface::makeShared<PerigeeSurface>(Vector3(0., 0., 0.));

    // Calculate d0 and z0 corresponding to vertex position
    double d0V = sqrt(x * x + y * y);
    double z0V = z;

    // Start constructing nTracks tracks in the following
    std::vector<InputTrack> tracks;

    // Construct random track emerging from vicinity of vertex position
    // Vector to store track objects used for vertex fit
    for (unsigned int iTrack = 0; iTrack < nTracks; iTrack++) {
      // Construct positive or negative charge randomly
      double q = qDist(gen) < 0 ? -1. : 1.;

      // Construct random track parameters
      BoundVector paramVec;
      paramVec << d0V + d0Dist(gen), z0V + z0Dist(gen), phiDist(gen),
          thetaDist(gen), q / pTDist(gen), 0.;

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
      tracks.push_back(InputTrack(
          BoundTrackParameters(perigeeSurface, paramVec, std::move(covMat))));
    }

    std::vector<const InputTrack*> tracksPtr;
    for (const auto& trk : tracks) {
      tracksPtr.push_back(&trk);
    }

    // Do the actual fit with 4 tracks without constraint
    Vertex<InputTrack> fittedVertex =
        billoirFitter.fit(tracksPtr, linearizer, vfOptions, state).value();
    if (!fittedVertex.tracks().empty()) {
      CHECK_CLOSE_ABS(fittedVertex.position(), vertexPosition, 1_mm);
    }
    // Do the fit with a constraint
    Vertex<InputTrack> fittedVertexConstraint =
        billoirFitter.fit(tracksPtr, linearizer, vfOptionsConstr, state)
            .value();
    if (!fittedVertexConstraint.tracks().empty()) {
      CHECK_CLOSE_ABS(fittedVertexConstraint.position(), vertexPosition, 1_mm);
    }

    std::cout << "Fitting nTracks: " << nTracks << std::endl;
    std::cout << "True Vertex: " << x << ", " << y << ", " << z << std::endl;
    std::cout << "Fitted Vertex: " << fittedVertex.position() << std::endl;
    std::cout << "Fitted constraint Vertex: "
              << fittedVertexConstraint.position() << std::endl;
  }
}

}  // namespace Test
}  // namespace Acts

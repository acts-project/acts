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
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <random>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

using Covariance = BoundSquareMatrix;
using Linearizer = HelicalTrackLinearizer<Propagator<EigenStepper<>>>;

// Create a test context
GeometryContext geoContext = GeometryContext();
MagneticFieldContext magFieldContext = MagneticFieldContext();

// 4D vertex distributions
// x-/y-position
std::uniform_real_distribution<> vXYDist(-0.1_mm, 0.1_mm);
// z-position
std::uniform_real_distribution<> vZDist(-20_mm, 20_mm);
// time
std::uniform_real_distribution<> vTDist(-1_ns, 1_ns);

// Track parameter distributions
// d0
std::uniform_real_distribution<> d0Dist(-0.01_mm, 0.01_mm);
// z0
std::uniform_real_distribution<> z0Dist(-0.2_mm, 0.2_mm);
// pT
std::uniform_real_distribution<> pTDist(0.4_GeV, 10_GeV);
// phi
std::uniform_real_distribution<> phiDist(-M_PI, M_PI);
// theta
std::uniform_real_distribution<> thetaDist(1.0, M_PI - 1.0);
// charge helper
std::uniform_real_distribution<> qDist(-1, 1);
// time
std::uniform_real_distribution<> tDist(-0.002_ns, 0.002_ns);

// Track parameter resolution distributions
// impact parameters
std::uniform_real_distribution<> resIPDist(0., 100_um);
// angles
std::uniform_real_distribution<> resAngDist(0., 0.1);
// q/p
std::uniform_real_distribution<> resQoPDist(-0.1, 0.1);
// Track time resolution distribution
std::uniform_real_distribution<> resTDist(0.1_ns, 0.2_ns);

// Number of tracks distritbution
std::uniform_int_distribution<> nTracksDist(3, 10);

// Dummy user-defined InputTrack type
struct InputTrack {
  InputTrack(const BoundTrackParameters& params) : m_parameters(params) {}

  const BoundTrackParameters& parameters() const { return m_parameters; }

  // store e.g. link to original objects here

 private:
  BoundTrackParameters m_parameters;
};

///
/// @brief Unit test for FullBilloirVertexFitter
/// with default input track type (= BoundTrackParameters)
///
BOOST_AUTO_TEST_CASE(billoir_vertex_fitter_defaulttrack_test) {
  // Set up RNG
  int seed = 31415;
  std::mt19937 gen(seed);
  // Set up constant B-Field
  auto bField = std::make_shared<ConstantBField>(Vector3{0.0, 0.0, 1_T});

  // Set up Eigenstepper
  EigenStepper<> stepper(bField);
  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator<EigenStepper<>>>(stepper);

  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  // Constraint for vertex fit
  Vertex<BoundTrackParameters> constraint;
  Vertex<InputTrack> customConstraint;
  // Some arbitrary values
  SquareMatrix4 covMatVtx = SquareMatrix4::Zero();
  double ns2 = Acts::UnitConstants::ns * Acts::UnitConstants::ns;
  covMatVtx(0, 0) = 30_mm2;
  covMatVtx(1, 1) = 30_mm2;
  covMatVtx(2, 2) = 30_mm2;
  covMatVtx(3, 3) = 30 * ns2;
  constraint.setFullCovariance(covMatVtx);
  constraint.setFullPosition(Vector4(0, 0, 0, 0));
  customConstraint.setFullCovariance(covMatVtx);
  customConstraint.setFullPosition(Vector4(0, 0, 0, 0));

  // Set up Billoir vertex fitter with default tracks
  using VertexFitter =
      FullBilloirVertexFitter<BoundTrackParameters, Linearizer>;
  VertexFitter::Config vertexFitterCfg;
  VertexFitter billoirFitter(vertexFitterCfg);
  VertexFitter::State state(bField->makeCache(magFieldContext));
  // Vertexing options for default tracks
  VertexingOptions<BoundTrackParameters> vfOptions(geoContext, magFieldContext);
  VertexingOptions<BoundTrackParameters> vfOptionsConstr(
      geoContext, magFieldContext, constraint);

  // Create a custom std::function to extract BoundTrackParameters from
  // user-defined InputTrack
  std::function<BoundTrackParameters(InputTrack)> extractParameters =
      [](const InputTrack& params) { return params.parameters(); };

  // Set up Billoir vertex fitter with user-defined input tracks
  using CustomVertexFitter = FullBilloirVertexFitter<InputTrack, Linearizer>;
  CustomVertexFitter::Config customVertexFitterCfg;
  CustomVertexFitter customBilloirFitter(customVertexFitterCfg,
                                         extractParameters);
  CustomVertexFitter::State customState(bField->makeCache(magFieldContext));
  // Vertexing options for custom tracks
  VertexingOptions<InputTrack> customVfOptions(geoContext, magFieldContext);
  VertexingOptions<InputTrack> customVfOptionsConstr(
      geoContext, magFieldContext, customConstraint);

  BOOST_TEST_CONTEXT(
      "Testing FullBilloirVertexFitter when input track vector is empty.") {
    const std::vector<const BoundTrackParameters*> emptyVector;

    // Without constraint
    Vertex<BoundTrackParameters> fittedVertex =
        billoirFitter.fit(emptyVector, linearizer, vfOptions, state).value();

    Vector3 origin(0., 0., 0.);
    SquareMatrix4 zeroMat = SquareMatrix4::Zero();
    BOOST_CHECK_EQUAL(fittedVertex.position(), origin);
    BOOST_CHECK_EQUAL(fittedVertex.fullCovariance(), zeroMat);

    // With constraint
    fittedVertex =
        billoirFitter.fit(emptyVector, linearizer, vfOptionsConstr, state)
            .value();

    BOOST_CHECK_EQUAL(fittedVertex.position(), origin);
    BOOST_CHECK_EQUAL(fittedVertex.fullCovariance(), zeroMat);
  }

  // Number of events
  const int nEvents = 100;
  for (int eventIdx = 0; eventIdx < nEvents; ++eventIdx) {
    // Number of tracks
    unsigned int nTracks = nTracksDist(gen);

    // Create position of vertex and perigee surface
    double x = vXYDist(gen);
    double y = vXYDist(gen);
    double z = vZDist(gen);
    double t = vTDist(gen);

    Vector4 trueVertex(x, y, z, t);
    std::shared_ptr<PerigeeSurface> perigeeSurface =
        Surface::makeShared<PerigeeSurface>(Vector3(0., 0., 0.));

    // Calculate d0 and z0 corresponding to the vertex position
    double d0V = std::hypot(x, y);
    double z0V = z;

    // Vector to store track objects used for vertex fit
    std::vector<BoundTrackParameters> tracks;
    std::vector<InputTrack> customTracks;

    // Calculate random track emerging from vicinity of vertex position
    for (unsigned int iTrack = 0; iTrack < nTracks; iTrack++) {
      // Charge
      double q = qDist(gen) < 0 ? -1. : 1.;

      // Track parameters
      BoundVector paramVec;
      paramVec << d0V + d0Dist(gen), z0V + z0Dist(gen), phiDist(gen),
          thetaDist(gen), q / pTDist(gen), t + tDist(gen);

      // Resolutions
      double resD0 = resIPDist(gen);
      double resZ0 = resIPDist(gen);
      double resPh = resAngDist(gen);
      double resTh = resAngDist(gen);
      double resQp = resQoPDist(gen);
      double resT = resTDist(gen);

      // Fill vector of track objects with simple covariance matrix
      Covariance covMat;

      covMat << resD0 * resD0, 0., 0., 0., 0., 0., 0., resZ0 * resZ0, 0., 0.,
          0., 0., 0., 0., resPh * resPh, 0., 0., 0., 0., 0., 0., resTh * resTh,
          0., 0., 0., 0., 0., 0., resQp * resQp, 0., 0., 0., 0., 0., 0.,
          resT * resT;
      tracks.emplace_back(BoundTrackParameters(perigeeSurface, paramVec, covMat,
                                               ParticleHypothesis::pion()));
      customTracks.emplace_back(
          BoundTrackParameters(perigeeSurface, paramVec, std::move(covMat),
                               ParticleHypothesis::pion()));
    }

    std::vector<const BoundTrackParameters*> tracksPtr;
    for (const auto& trk : tracks) {
      tracksPtr.push_back(&trk);
    }

    std::vector<const InputTrack*> customTracksPtr;
    for (const auto& trk : customTracks) {
      customTracksPtr.push_back(&trk);
    }

    auto fit = [&trueVertex, &nTracks](const auto& fitter, const auto& trksPtr,
                                       const auto& lin, const auto& vfOpts,
                                       auto& vfState) {
      auto fittedVertex = fitter.fit(trksPtr, lin, vfOpts, vfState).value();
      if (!fittedVertex.tracks().empty()) {
        CHECK_CLOSE_ABS(fittedVertex.position(), trueVertex.head(3), 1_mm);
        auto tracksAtVtx = fittedVertex.tracks();
        auto trackAtVtx = tracksAtVtx[0];
        CHECK_CLOSE_ABS(fittedVertex.time(), trueVertex[3], 1_ns);
      }

      std::cout << "\nFitting " << nTracks << " tracks" << std::endl;
      std::cout << "True Vertex:\n" << trueVertex << std::endl;
      std::cout << "Fitted Vertex:\n"
                << fittedVertex.fullPosition() << std::endl;
    };

    BOOST_TEST_CONTEXT(
        "Testing FullBilloirVertexFitter without vertex constraint.") {
      fit(billoirFitter, tracksPtr, linearizer, vfOptions, state);
    }
    BOOST_TEST_CONTEXT(
        "Testing FullBilloirVertexFitter with vertex constraint.") {
      fit(billoirFitter, tracksPtr, linearizer, vfOptionsConstr, state);
    }
    BOOST_TEST_CONTEXT(
        "Testing FullBilloirVertexFitter with custom tracks (no vertex "
        "constraint).") {
      fit(customBilloirFitter, customTracksPtr, linearizer, customVfOptions,
          customState);
    }
    BOOST_TEST_CONTEXT(
        "Testing FullBilloirVertexFitter with custom tracks (with vertex "
        "constraint).") {
      fit(customBilloirFitter, customTracksPtr, linearizer,
          customVfOptionsConstr, customState);
    }
  }
}
}  // namespace Test
}  // namespace Acts

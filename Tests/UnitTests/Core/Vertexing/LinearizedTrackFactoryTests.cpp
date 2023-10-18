// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
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
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"
#include "Acts/Vertexing/NumericalTrackLinearizer.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <random>
#include <tuple>
#include <utility>
#include <vector>

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

using Covariance = BoundSquareMatrix;
// We will compare analytical and numerical computations in the case of a
// (non-zero) constant B-field and a zero B-field.
using HelicalPropagator = Propagator<EigenStepper<>>;
using StraightPropagator = Propagator<StraightLineStepper>;
using AnalyticalLinearizer = HelicalTrackLinearizer<HelicalPropagator>;
using StraightAnalyticalLinearizer = HelicalTrackLinearizer<StraightPropagator>;
using NumericalLinearizer = NumericalTrackLinearizer<HelicalPropagator>;
using StraightNumericalLinearizer =
    NumericalTrackLinearizer<StraightPropagator>;

// Create a test context
GeometryContext geoContext = GeometryContext();
MagneticFieldContext magFieldContext = MagneticFieldContext();

// Vertex x/y position distribution
std::uniform_real_distribution<> vXYDist(-0.1_mm, 0.1_mm);
// Vertex z position distribution
std::uniform_real_distribution<> vZDist(-20_mm, 20_mm);
// Vertex time distribution
std::uniform_real_distribution<> vTDist(-1_ns, 1_ns);
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
// Track time distribution
std::uniform_real_distribution<> tDist(-0.002_ns, 0.002_ns);
// Track IP resolution distribution
std::uniform_real_distribution<> resIPDist(0., 100_um);
// Track angular distribution
std::uniform_real_distribution<> resAngDist(0., 0.1);
// Track q/p resolution distribution
std::uniform_real_distribution<> resQoPDist(0.0, 0.1);
// Track time resolution distribution
std::uniform_real_distribution<> resTDist(0.1_ns, 0.2_ns);

///
/// @brief Test HelicalTrackLinearizer by comparing it to NumericalTrackLinearizer.
///
/// @note While HelicalTrackLinearizer computes the Jacobians using analytically derived formulae,
/// NumericalTrackLinearizer uses numerical differentiation:
/// f'(x) ~ (f(x+dx) - f(x)) / dx).
///
BOOST_AUTO_TEST_CASE(linearized_track_factory_test) {
  // Number of tracks to linearize
  unsigned int nTracks = 100;

  // Set up RNG
  int seed = 31415;
  std::mt19937 gen(seed);

  // Constant B-Field and 0 B-field
  auto constField = std::make_shared<ConstantBField>(Vector3{0.0, 0.0, 2_T});
  auto zeroField = std::make_shared<NullBField>();

  // Set up stepper and propagator for constant B-field
  EigenStepper<> stepper(constField);
  auto propagator = std::make_shared<HelicalPropagator>(stepper);

  // Set up stepper and propagator for 0 B-field
  StraightLineStepper straightStepper;
  auto straightPropagator =
      std::make_shared<StraightPropagator>(straightStepper);

  // Create perigee surface, initial track parameters will be relative to it
  std::shared_ptr<PerigeeSurface> perigeeSurface{
      Surface::makeShared<PerigeeSurface>(Vector3{0., 0., 0.})};

  // Vertex position and corresponding d0 and z0
  Vector4 vtxPos;
  double d0v{};
  double z0v{};
  double t0v{};
  {
    double x = vXYDist(gen);
    double y = vXYDist(gen);
    double z = vZDist(gen);
    double t = vTDist(gen);
    d0v = std::hypot(x, y);
    z0v = z;
    t0v = t;
    vtxPos << x, y, z, t;
  }

  // Vector storing the tracks that we linearize
  std::vector<BoundTrackParameters> tracks;

  // Construct random track emerging from vicinity of vertex position
  for (unsigned int iTrack = 0; iTrack < nTracks; iTrack++) {
    // Random charge
    double q = qDist(gen) < 0 ? -1. : 1.;

    // Random track parameters
    BoundVector paramVec;
    paramVec << d0v + d0Dist(gen), z0v + z0Dist(gen), phiDist(gen),
        thetaDist(gen), q / pTDist(gen), t0v + tDist(gen);

    // Resolutions
    double resD0 = resIPDist(gen);
    double resZ0 = resIPDist(gen);
    double resPh = resAngDist(gen);
    double resTh = resAngDist(gen);
    double resQp = resQoPDist(gen);
    double resT = resTDist(gen);

    // Fill vector of track objects with simple covariance matrix
    Covariance covMat;

    covMat << resD0 * resD0, 0., 0., 0., 0., 0., 0., resZ0 * resZ0, 0., 0., 0.,
        0., 0., 0., resPh * resPh, 0., 0., 0., 0., 0., 0., resTh * resTh, 0.,
        0., 0., 0., 0., 0., resQp * resQp, 0., 0., 0., 0., 0., 0., resT * resT;
    tracks.emplace_back(perigeeSurface, paramVec, std::move(covMat),
                        ParticleHypothesis::pion());
  }

  // Linearizer for constant field and corresponding state
  AnalyticalLinearizer::Config linConfig(constField, propagator);
  AnalyticalLinearizer linFactory(linConfig);
  AnalyticalLinearizer::State linState(constField->makeCache(magFieldContext));

  NumericalLinearizer::Config numLinConfig(constField, propagator);
  NumericalLinearizer numLinFactory(numLinConfig);
  NumericalLinearizer::State numLinState(
      constField->makeCache(magFieldContext));

  // Linearizer for 0 field and corresponding state
  StraightAnalyticalLinearizer::Config straightLinConfig(straightPropagator);
  StraightAnalyticalLinearizer straightLinFactory(straightLinConfig);
  StraightAnalyticalLinearizer::State straightLinState(
      zeroField->makeCache(magFieldContext));

  StraightNumericalLinearizer::Config numStraightLinConfig(straightPropagator);
  StraightNumericalLinearizer numStraightLinFactory(numStraightLinConfig);
  StraightNumericalLinearizer::State numStraightLinState(
      zeroField->makeCache(magFieldContext));

  // Lambda for comparing outputs of the two linearization methods
  // We compare the linearization result at the PCA to "linPoint"
  auto checkLinearizers = [](auto& lin1, auto& linState1, auto& lin2,
                             auto& linState2, const BoundTrackParameters& track,
                             const Vector4& linPoint,
                             const auto& geometryContext,
                             const auto& fieldContext) {
    // In addition to comparing the output of the linearizers, we check that
    // they return non-zero quantities
    BoundVector vecBoundZero = BoundVector::Zero();
    BoundSquareMatrix matBoundZero = BoundSquareMatrix::Zero();
    ActsMatrix<eBoundSize, 4> matBound2SPZero =
        ActsMatrix<eBoundSize, 4>::Zero();
    ActsMatrix<eBoundSize, 3> matBound2MomZero =
        ActsMatrix<eBoundSize, 3>::Zero();

    // We check that the entries of the output quantities either
    // -) have a relative difference of less than "relTol"
    // or
    // -) are both smaller than "small"
    double relTol = 5e-4;
    double small = 5e-4;

    std::shared_ptr<PerigeeSurface> perigee =
        Surface::makeShared<PerigeeSurface>(VectorHelpers::position(linPoint));

    const LinearizedTrack linTrack1 =
        lin1.linearizeTrack(track, linPoint[3], *perigee, geometryContext,
                            fieldContext, linState1)
            .value();
    const LinearizedTrack linTrack2 =
        lin2.linearizeTrack(track, linPoint[3], *perigee, geometryContext,
                            fieldContext, linState2)
            .value();

    // There should be no problem here because both linearizers compute
    // "parametersAtPCA" the same way
    CHECK_CLOSE_OR_SMALL(linTrack1.parametersAtPCA, linTrack2.parametersAtPCA,
                         relTol, small);
    BOOST_CHECK_NE(linTrack1.parametersAtPCA, vecBoundZero);
    BOOST_CHECK_NE(linTrack2.parametersAtPCA, vecBoundZero);

    // Compare position Jacobians
    CHECK_CLOSE_OR_SMALL((linTrack1.positionJacobian),
                         (linTrack2.positionJacobian), relTol, small);
    BOOST_CHECK_NE(linTrack1.positionJacobian, matBound2SPZero);
    BOOST_CHECK_NE(linTrack2.positionJacobian, matBound2SPZero);

    // Compare momentum Jacobians
    CHECK_CLOSE_OR_SMALL((linTrack1.momentumJacobian),
                         (linTrack2.momentumJacobian), relTol, small);
    BOOST_CHECK_NE(linTrack1.momentumJacobian, matBound2MomZero);
    BOOST_CHECK_NE(linTrack2.momentumJacobian, matBound2MomZero);

    // Again, both methods compute "covarianceAtPCA" the same way => this
    // check should always work
    CHECK_CLOSE_OR_SMALL(linTrack1.covarianceAtPCA, linTrack2.covarianceAtPCA,
                         relTol, small);
    BOOST_CHECK_NE(linTrack1.covarianceAtPCA, matBoundZero);
    BOOST_CHECK_NE(linTrack2.covarianceAtPCA, matBoundZero);

    // Check whether "linPoint" is saved correctly in the LinearizerTrack
    // objects
    BOOST_CHECK_EQUAL(linTrack1.linearizationPoint, linPoint);
    BOOST_CHECK_EQUAL(linTrack2.linearizationPoint, linPoint);

    CHECK_CLOSE_OR_SMALL(linTrack1.constantTerm, linTrack2.constantTerm, relTol,
                         small);
    BOOST_CHECK_NE(linTrack1.constantTerm, vecBoundZero);
    BOOST_CHECK_NE(linTrack2.constantTerm, vecBoundZero);
  };

  // Compare linearizers for all tracks
  for (const BoundTrackParameters& trk : tracks) {
    BOOST_TEST_CONTEXT("Linearization in constant magnetic field") {
      checkLinearizers(linFactory, linState, numLinFactory, numLinState, trk,
                       vtxPos, geoContext, magFieldContext);
    }
    BOOST_TEST_CONTEXT("Linearization without magnetic field") {
      checkLinearizers(straightLinFactory, straightLinState,
                       numStraightLinFactory, numStraightLinState, trk, vtxPos,
                       geoContext, magFieldContext);
    }
  }
}

}  // namespace Test
}  // namespace Acts

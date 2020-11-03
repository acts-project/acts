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

#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/Vertex.hpp"

#include <memory>
#include <utility>

namespace {

using namespace Acts;
using namespace Acts::UnitLiterals;
using Acts::VectorHelpers::makeVector4;

using MagneticField = Acts::ConstantBField;
using Stepper = Acts::EigenStepper<MagneticField>;
using Propagator = Acts::Propagator<Stepper>;

using Estimator =
    Acts::ImpactPointEstimator<Acts::BoundTrackParameters, Propagator>;

const Acts::GeometryContext geoContext;
const Acts::MagneticFieldContext magFieldContext;

// use a fixed seed for reproducible tests
std::default_random_engine gen(31415);

// Construct an impact point estimator for a constant bfield along z.
Estimator makeEstimator(double bZ) {
  MagneticField field(Vector3D(0, 0, bZ));
  Stepper stepper(field);
  Estimator::Config cfg(field,
                        std::make_shared<Propagator>(std::move(stepper)));
  return Estimator(cfg);
}

// Generate track parameter vector and diagonal covariance matrix.
template <typename rng_t>
std::pair<Acts::BoundVector, Acts::BoundSymMatrix>
makeTrackParametersCovariance(rng_t& rng) {
  using namespace Acts;
  using Uniform = std::uniform_real_distribution<BoundScalar>;

  auto boolDist = std::uniform_int_distribution<int>(0, 1);

  // generate parameters uniformly within reasonable ranges
  BoundVector params;
  params[eBoundLoc0] = Uniform(-10_um, 10_um)(rng);
  params[eBoundLoc1] = Uniform(-200_um, 200_um)(rng);
  params[eBoundTime] = Uniform(-5_ns, 5_ns)(rng);
  params[eBoundPhi] = Uniform(-M_PI, M_PI)(rng);
  params[eBoundTheta] = Uniform(1.0, M_PI - 1.0)(rng);
  params[eBoundQOverP] =
      (boolDist(rng) ? -1_e : 1_e) / Uniform(0.4_GeV, 10_GeV)(rng);
  // generate random resolutions w/p correlations in reasonable ranges
  // this ignores correlations between parameter values and their resolution
  BoundVector stddev;
  stddev[eBoundLoc0] = Uniform(5_um, 100_um)(rng);
  stddev[eBoundLoc1] = Uniform(5_um, 100_um)(rng);
  stddev[eBoundTime] = Uniform(1_ns, 5_ns)(rng);
  stddev[eBoundPhi] = Uniform(0.1_degree, 1_degree)(rng);
  stddev[eBoundTheta] = Uniform(0.1_degree, 1_degree)(rng);
  // draw relative momentum resolution
  stddev[eBoundQOverP] = Uniform(0.0125, 0.125)(rng) * params[eBoundQOverP];

  return {params, stddev.cwiseProduct(stddev).asDiagonal()};
}

// Generate vertex position vector and diagonal covariance matrix.
template <typename rng_t>
std::pair<Acts::Vector4D, Acts::SymMatrix4D> makeVertexParametersCovariance(
    rng_t& rng) {
  using namespace Acts;
  using Normal = std::normal_distribution<double>;
  using Uniform = std::uniform_real_distribution<BoundScalar>;

  auto stdNormalDist = std::normal_distribution<double>(0, 1);

  Vector4D stddev;
  stddev[ePos0] = Uniform(5_um, 50_um)(rng);
  stddev[ePos1] = Uniform(5_um, 50_um)(rng);
  stddev[ePos2] = Uniform(10_um, 100_um)(rng);
  stddev[eTime] = Uniform(1_ns, 5_ns)(rng);
  Vector4D pos;
  pos[ePos0] = stdNormalDist(rng) * stddev[ePos0];
  pos[ePos1] = stdNormalDist(rng) * stddev[ePos1];
  pos[ePos2] = stdNormalDist(rng) * stddev[ePos2] + Uniform(-20_mm, 20_mm)(rng);
  pos[eTime] = stdNormalDist(rng) * stddev[eTime];

  return {pos, stddev.cwiseProduct(stddev).asDiagonal()};
}

}  // namespace

BOOST_AUTO_TEST_SUITE(Vertexing)

/// @brief Unit test for ImpactPointEstimator params and distance
///
BOOST_AUTO_TEST_CASE(ImpactPointEstimator3d) {
  // Debug mode
  bool debugMode = false;
  // Number of tests
  unsigned int nTests = 10;

  Estimator ipEstimator = makeEstimator(2_T);
  Estimator::State state(magFieldContext);

  // Reference position and corresponding perigee surface
  Vector3D refPosition(0., 0., 0.);
  auto perigeeSurface = Surface::makeShared<PerigeeSurface>(refPosition);

  // Start running tests
  for (unsigned int i = 0; i < nTests; i++) {
    // Creating the track
    auto [par, cov] = makeTrackParametersCovariance(gen);
    BoundTrackParameters myTrack(perigeeSurface, par, cov);

    // Estimate 3D distance must be less than the 2d distance on the perigee
    double distT = std::hypot(par[eBoundLoc0], par[eBoundLoc1]);
    double dist3 =
        ipEstimator.calculate3dDistance(geoContext, myTrack, refPosition, state)
            .value();
    BOOST_CHECK_LT(dist3, distT);

    if (debugMode) {
      std::cout << std::setprecision(10)
                << "Distance in transverse plane: " << distT << std::endl;
      std::cout << std::setprecision(10) << "Distance in 3D: " << dist3
                << std::endl;
    }

    auto res = ipEstimator.estimate3DImpactParameters(
        geoContext, magFieldContext, myTrack, refPosition, state);
    BoundTrackParameters trackAtIP3d = **res;
    const auto& myTrackParams = myTrack.parameters();
    const auto& trackIP3dParams = trackAtIP3d.parameters();

    // d0 and z0 should have changed
    BOOST_CHECK_NE(myTrackParams[eBoundLoc0], trackIP3dParams[eBoundLoc0]);
    BOOST_CHECK_NE(myTrackParams[eBoundLoc1], trackIP3dParams[eBoundLoc1]);
    // Theta along helix and q/p shoud remain the same
    CHECK_CLOSE_REL(myTrackParams[eBoundTheta], trackIP3dParams[eBoundTheta],
                    1e-5);
    CHECK_CLOSE_REL(myTrackParams[eBoundQOverP], trackIP3dParams[eBoundQOverP],
                    1e-5);

    if (debugMode) {
      std::cout << std::setprecision(10) << "Old track parameters: \n"
                << myTrackParams << std::endl;
      std::cout << std::setprecision(10) << "Parameters at IP3d: \n"
                << trackIP3dParams << std::endl;
    }
  }  // end for loop tests
}

/// @brief Unit test for ImpactPointEstimator
///  compatibility estimator
BOOST_AUTO_TEST_CASE(impactpoint_estimator_compatibility_test) {
  // Debug mode
  bool debugMode = false;
  // Number of tests
  unsigned int nTests = 10;

  Estimator ipEstimator = makeEstimator(2_T);
  Estimator::State state(magFieldContext);

  // Reference position and corresponding perigee surface
  Vector3D refPosition(0., 0., 0.);
  auto perigeeSurface = Surface::makeShared<PerigeeSurface>(refPosition);

  // Lists to store distances and comp values
  std::vector<double> distancesList;
  std::vector<double> compatibilityList;

  // Start running tests
  for (unsigned int i = 0; i < nTests; i++) {
    // Creating the track
    auto [par, cov] = makeTrackParametersCovariance(gen);
    BoundTrackParameters myTrack(perigeeSurface, par, cov);

    // Estimate 3D distance
    auto distanceRes = ipEstimator.calculate3dDistance(geoContext, myTrack,
                                                       refPosition, state);
    BOOST_CHECK(distanceRes.ok());

    distancesList.push_back(*distanceRes);

    auto res = ipEstimator.estimate3DImpactParameters(
        geoContext, magFieldContext, myTrack, refPosition, state);

    BOOST_CHECK(res.ok());

    BoundTrackParameters params = std::move(**res);

    auto compRes =
        ipEstimator.get3dVertexCompatibility(geoContext, &params, refPosition);

    BOOST_CHECK(compRes.ok());

    compatibilityList.push_back(*compRes);

  }  // end create tracks loop

  // Now test for all above constructed tracks
  // if distances and compatibility values are
  // compatible with one another
  for (unsigned int i = 0; i < nTests; i++) {
    for (unsigned int j = i + 1; j < nTests; j++) {
      double relDiffComp =
          (compatibilityList[i] - compatibilityList[j]) / compatibilityList[i];

      double relDiffDist =
          (distancesList[i] - distancesList[j]) / distancesList[i];

      if (debugMode) {
        std::cout << "Comparing track " << i << " with track " << j
                  << std::endl;
        std::cout << "\t" << i << ": Comp.: " << compatibilityList[i]
                  << ", dist.: " << distancesList[i] << std::endl;
        std::cout << "\t" << j << ": Comp.: " << compatibilityList[j]
                  << ", dist.: " << distancesList[j] << std::endl;
        std::cout << "\t Rel.diff.: Comp(1-2)/1: " << relDiffComp
                  << ", Dist(1-2)/1: " << relDiffDist << std::endl;
      }

      // Relative differences of compatibility values and distances
      // should have the same sign, i.e. the following product
      // should always be positive
      double res = relDiffComp * relDiffDist;

      // TODO 2020-09-09 msmk
      // this fails for one track after the the track parameters cleanup.
      // i do not understand what this tests and/or how to fix it. Bastian
      // has to look at this.
      BOOST_CHECK_GE(res, 0);
    }
  }
}

/// @brief Unit test for ImpactPoint 3d estimator, using same
/// configuration and test values as in Athena unit test algorithm
/// Tracking/TrkVertexFitter/TrkVertexFitterUtils/test/ImpactPointEstimator_test
BOOST_AUTO_TEST_CASE(impactpoint_estimator_athena_test) {
  Estimator ipEstimator = makeEstimator(1.9971546939_T);
  Estimator::State state(magFieldContext);

  // Use same values as in Athena unit test
  Vector4D pos1(2_mm, 1_mm, -10_mm, 0_ns);
  Vector3D mom1(400_MeV, 600_MeV, 200_MeV);
  Vector3D vtxPos(1.2_mm, 0.8_mm, -7_mm);

  // Start creating some track parameters
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(pos1.segment<3>(ePos0));
  // Some fixed track parameter values
  BoundTrackParameters params1(
      perigeeSurface, geoContext, pos1, mom1, mom1.norm(), 1_e,
      BoundTrackParameters::CovarianceMatrix::Identity());

  auto res1 =
      ipEstimator.calculate3dDistance(geoContext, params1, vtxPos, state);
  BOOST_CHECK(res1.ok());
  double distance = (*res1);

  // Desired result from Athena unit test
  CHECK_CLOSE_ABS(distance, 3.10391_mm, 0.00001_mm);

  auto res2 = ipEstimator.estimate3DImpactParameters(
      geoContext, magFieldContext, params1, vtxPos, state);
  BOOST_CHECK(res2.ok());
  BoundTrackParameters endParams = std::move(**res2);
  Vector3D surfaceCenter = endParams.referenceSurface().center(geoContext);

  BOOST_CHECK_EQUAL(surfaceCenter, vtxPos);
}

///
/// @brief Unit test for impact parameter estimation
///
BOOST_AUTO_TEST_CASE(impactpoint_estimator_parameter_estimation_test) {
  // Number of tracks to test with
  unsigned int nTracks = 10;

  Estimator ipEstimator = makeEstimator(1_T);
  Estimator::State state(magFieldContext);

  // Reference position and corresponding perigee surface
  Vector3D refPosition(0., 0., 0.);
  auto perigeeSurface = Surface::makeShared<PerigeeSurface>(refPosition);

  // Create position of vertex and perigee surface
  auto [vtxPos, vtxCov] = makeVertexParametersCovariance(gen);
  Vertex<BoundTrackParameters> myConstraint;
  myConstraint.setFullPosition(vtxPos);
  myConstraint.setFullCovariance(4 * vtxCov);

  // Construct random track emerging from vicinity of vertex position
  // Vector to store track objects used for vertex fit
  for (unsigned int iTrack = 0; iTrack < nTracks; iTrack++) {
    auto [par, cov] = makeTrackParametersCovariance(gen);
    BoundTrackParameters track(perigeeSurface, par, cov);

    // Check if IP are retrieved
    ImpactParametersAndSigma output =
        ipEstimator
            .estimateImpactParameters(track, myConstraint, geoContext,
                                      magFieldContext)
            .value();
    BOOST_CHECK_NE(output.IPd0, 0.);
    BOOST_CHECK_NE(output.IPz0, 0.);
  }
}

BOOST_AUTO_TEST_SUITE_END()

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
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/ImpactPoint3dEstimator.hpp"
#include "Acts/Vertexing/Vertex.hpp"

using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

using Covariance = BoundSymMatrix;
using Propagator = Propagator<EigenStepper<ConstantBField>>;
// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();

// Track d0 distribution
std::uniform_real_distribution<> d0Dist(-0.01_mm, 0.01_mm);
// Track z0 distribution
std::uniform_real_distribution<> z0Dist(-0.2_mm, 0.2_mm);
// Track pT distribution
std::uniform_real_distribution<> pTDist(0.4_GeV, 10._GeV);
// Track phi distribution
std::uniform_real_distribution<> phiDist(-M_PI, M_PI);
// Track theta distribution
std::uniform_real_distribution<> thetaDist(1.0, M_PI - 1.0);
// Track IP resolution distribution
std::uniform_real_distribution<> resIPDist(0., 100._um);
// Track angular distribution
std::uniform_real_distribution<> resAngDist(0., 0.1);
// Track q/p resolution distribution
std::uniform_real_distribution<> resQoPDist(-0.1, 0.1);
// Track charge helper distribution
std::uniform_real_distribution<> qDist(-1, 1);

/// @brief Unit test for ImpactPoint3dEstimator params and distance
///
BOOST_AUTO_TEST_CASE(impactpoint_3d_estimator_params_distance_test) {
  // Debug mode
  bool debugMode = false;
  // Number of tests
  unsigned int nTests = 10;

  // Set up RNG
  int mySeed = 31415;
  std::mt19937 gen(mySeed);

  // Set up constant B-Field
  ConstantBField bField(Vector3D(0., 0., 2._T));

  // Set up Eigenstepper
  EigenStepper<ConstantBField> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);
  PropagatorOptions<> pOptions(tgContext, mfContext);

  // Set up the ImpactPoint3dEstimator
  ImpactPoint3dEstimator<BoundParameters, Propagator>::Config ipEstCfg(
      bField, propagator, pOptions);

  ImpactPoint3dEstimator<BoundParameters, Propagator> ipEstimator(ipEstCfg);

  // Reference position
  Vector3D refPosition(0., 0., 0.);

  // Start running tests
  for (unsigned int i = 0; i < nTests; i++) {
    // Create a track
    // Resolutions
    double resD0 = resIPDist(gen);
    double resZ0 = resIPDist(gen);
    double resPh = resAngDist(gen);
    double resTh = resAngDist(gen);
    double resQp = resQoPDist(gen);

    // Covariance matrix
    Covariance covMat;
    covMat << resD0 * resD0, 0., 0., 0., 0., 0., 0., resZ0 * resZ0, 0., 0., 0.,
        0., 0., 0., resPh * resPh, 0., 0., 0., 0., 0., 0., resTh * resTh, 0.,
        0., 0., 0., 0., 0., resQp * resQp, 0., 0., 0., 0., 0., 0., 1.;

    // The charge
    double q = qDist(gen) < 0 ? -1. : 1.;

    // Impact parameters (IP)
    double d0 = d0Dist(gen);
    double z0 = z0Dist(gen);

    if (debugMode) {
      std::cout << "IP: (" << d0 << "," << z0 << ")" << std::endl;
    }

    // The track parameters
    BoundParameters::ParVector_t paramVec;
    paramVec << d0, z0, phiDist(gen), thetaDist(gen), q / pTDist(gen), 0.;

    // Corresponding surface
    std::shared_ptr<PerigeeSurface> perigeeSurface =
        Surface::makeShared<PerigeeSurface>(refPosition);

    // Creating the track
    BoundParameters myTrack(tgContext, std::move(covMat), paramVec,
                            perigeeSurface);

    // Distance in transverse plane
    double transverseDist = std::sqrt(std::pow(d0, 2) + std::pow(z0, 2));

    // Estimate 3D distance
    auto distanceRes =
        ipEstimator.calculateDistance(tgContext, myTrack, refPosition);
    BOOST_CHECK(distanceRes.ok());

    double distance = *distanceRes;

    BOOST_CHECK(distance < transverseDist);

    if (debugMode) {
      std::cout << std::setprecision(10)
                << "Distance in transverse plane: " << transverseDist
                << std::endl;
      std::cout << std::setprecision(10) << "Distance in 3D: " << distance
                << std::endl;
    }

    auto res =
        ipEstimator.getParamsAtClosestApproach(tgContext, myTrack, refPosition);

    BOOST_CHECK(res.ok());

    BoundParameters trackAtIP3d = std::move(**res);

    const auto& myTrackParams = myTrack.parameters();
    const auto& trackIP3dParams = trackAtIP3d.parameters();

    // d0 and z0 should have changed
    BOOST_CHECK_NE(myTrackParams[ParID_t::eLOC_D0],
                   trackIP3dParams[ParID_t::eLOC_D0]);
    BOOST_CHECK_NE(myTrackParams[ParID_t::eLOC_Z0],
                   trackIP3dParams[ParID_t::eLOC_Z0]);
    // Theta along helix and q/p shoud remain the same
    CHECK_CLOSE_REL(myTrackParams[ParID_t::eTHETA],
                    trackIP3dParams[ParID_t::eTHETA], 1e-5);
    CHECK_CLOSE_REL(myTrackParams[ParID_t::eQOP],
                    trackIP3dParams[ParID_t::eQOP], 1e-5);

    if (debugMode) {
      std::cout << std::setprecision(10) << "Old track parameters: \n"
                << myTrackParams << std::endl;
      std::cout << std::setprecision(10) << "Parameters at IP3d: \n"
                << trackIP3dParams << std::endl;
    }
  }  // end for loop tests
}

/// @brief Unit test for ImpactPoint3dEstimator
///  compatibility estimator
BOOST_AUTO_TEST_CASE(impactpoint_3d_estimator_compatibility_test) {
  // Debug mode
  bool debugMode = false;
  // Number of tests
  unsigned int nTests = 10;

  // Set up RNG
  int mySeed = 31415;
  std::mt19937 gen(mySeed);

  // Set up constant B-Field
  ConstantBField bField(Vector3D(0., 0., 2.) * units::_T);

  // Set up Eigenstepper
  EigenStepper<ConstantBField> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);
  PropagatorOptions<> pOptions(tgContext, mfContext);

  // Set up the ImpactPoint3dEstimator
  ImpactPoint3dEstimator<BoundParameters, Propagator>::Config ipEstCfg(
      bField, propagator, pOptions);

  ImpactPoint3dEstimator<BoundParameters, Propagator> ipEstimator(ipEstCfg);

  // Reference position
  Vector3D refPosition(0., 0., 0.);

  // Lists to store distances and comp values
  std::vector<double> distancesList;
  std::vector<double> compatibilityList;

  // Generate covariance matrix values once
  // for all cov matrices below
  // Resolutions
  double resD0 = resIPDist(gen);
  double resZ0 = resIPDist(gen);
  double resPh = resAngDist(gen);
  double resTh = resAngDist(gen);
  double resQp = resQoPDist(gen);

  // Start running tests
  for (unsigned int i = 0; i < nTests; i++) {
    // Create a track

    // Covariance matrix
    Covariance covMat;
    covMat << resD0 * resD0, 0., 0., 0., 0., 0., 0., resZ0 * resZ0, 0., 0., 0.,
        0., 0., 0., resPh * resPh, 0., 0., 0., 0., 0., 0., resTh * resTh, 0.,
        0., 0., 0., 0., 0., resQp * resQp, 0., 0., 0., 0., 0., 0., 1.;

    // The charge
    double q = qDist(gen) < 0 ? -1. : 1.;

    // Impact parameters (IP)
    double d0 = d0Dist(gen);
    double z0 = z0Dist(gen);

    // The track parameters
    BoundParameters::ParVector_t paramVec;
    paramVec << d0, z0, phiDist(gen), thetaDist(gen), q / pTDist(gen), 0.;

    // Corresponding surface
    std::shared_ptr<PerigeeSurface> perigeeSurface =
        Surface::makeShared<PerigeeSurface>(refPosition);

    // Creating the track
    BoundParameters myTrack(tgContext, std::move(covMat), paramVec,
                            perigeeSurface);

    // Estimate 3D distance
    auto distanceRes =
        ipEstimator.calculateDistance(tgContext, myTrack, refPosition);
    BOOST_CHECK(distanceRes.ok());

    distancesList.push_back(*distanceRes);

    auto res =
        ipEstimator.getParamsAtClosestApproach(tgContext, myTrack, refPosition);

    BOOST_CHECK(res.ok());

    BoundParameters params = std::move(**res);

    auto compRes =
        ipEstimator.getVertexCompatibility(tgContext, &params, refPosition);

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

      BOOST_CHECK(res > 0.);
    }
  }
}

/// @brief Unit test for ImpactPoint 3d estimator, using same
/// configuration and test values as in Athena unit test algorithm
/// Tracking/TrkVertexFitter/TrkVertexFitterUtils/test/ImpactPoint3dEstimator_test
BOOST_AUTO_TEST_CASE(impactpoint_3d_estimator_athena_test) {
  // Set up constant B-Field
  ConstantBField bField(Vector3D(0., 0., 1.9971546939_T));

  // Set up Eigenstepper
  EigenStepper<ConstantBField> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);
  PropagatorOptions<> pOptions(tgContext, mfContext);

  // Set up the ImpactPoint3dEstimator
  ImpactPoint3dEstimator<BoundParameters, Propagator>::Config ipEstCfg(
      bField, propagator, pOptions);

  ImpactPoint3dEstimator<BoundParameters, Propagator> ipEstimator(ipEstCfg);

  // Use same values as in Athena unit test
  Vector3D pos1(2_mm, 1_mm, -10_mm);
  Vector3D mom1(400_MeV, 600_MeV, 200_MeV);
  Vector3D vtxPos(1.2_mm, 0.8_mm, -7_mm);

  // Start creating some track parameters
  Covariance covMat = Covariance::Identity();
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(pos1);

  // Some fixed track parameter values
  BoundParameters params1(tgContext, covMat, pos1, mom1, 1, 0, perigeeSurface);

  auto res1 = ipEstimator.calculateDistance(tgContext, params1, vtxPos);
  BOOST_CHECK(res1.ok());
  double distance = (*res1);

  // Desired result from Athena unit test
  const double result = 3.10391_mm;
  CHECK_CLOSE_ABS(distance, result, 0.00001_mm);

  auto res2 =
      ipEstimator.getParamsAtClosestApproach(tgContext, params1, vtxPos);
  BOOST_CHECK(res2.ok());
  BoundParameters endParams = std::move(**res2);
  Vector3D surfaceCenter = endParams.referenceSurface().center(tgContext);

  BOOST_CHECK_EQUAL(surfaceCenter, vtxPos);
}

}  // namespace Test
}  // namespace Acts

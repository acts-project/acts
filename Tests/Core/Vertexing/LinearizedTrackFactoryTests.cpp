// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE LinearizedTrackFactory Tests
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/LinearizedTrackFactory.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"

namespace bdata = boost::unit_test::data;

namespace Acts {
namespace Test {

  // Vertex x/y position distribution
  std::uniform_real_distribution<> vXYDist(-0.1 * units::_mm, 0.1 * units::_mm);
  // Vertex z position distribution
  std::uniform_real_distribution<> vZDist(-20 * units::_mm, 20 * units::_mm);
  // Track d0 distribution
  std::uniform_real_distribution<> d0Dist(-0.01 * units::_mm,
                                          0.01 * units::_mm);
  // Track z0 distribution
  std::uniform_real_distribution<> z0Dist(-0.2 * units::_mm, 0.2 * units::_mm);
  // Track pT distribution
  std::uniform_real_distribution<> pTDist(0.4 * units::_GeV, 10. * units::_GeV);
  // Track phi distribution
  std::uniform_real_distribution<> phiDist(-M_PI, M_PI);
  // Track theta distribution
  std::uniform_real_distribution<> thetaDist(1.0, M_PI - 1.0);
  // Track charge helper distribution
  std::uniform_real_distribution<> qDist(-1, 1);
  // Track IP resolution distribution
  std::uniform_real_distribution<> resIPDist(0., 100. * units::_um);
  // Track angular distribution
  std::uniform_real_distribution<> resAngDist(0., 0.1);
  // Track q/p resolution distribution
  std::uniform_real_distribution<> resQoPDist(-0.1, 0.1);

  ///
  /// @brief Unit test for LinearizedTrackFactory
  ///
  BOOST_AUTO_TEST_CASE(linearized_track_factory_test)
  {

    // Number of tracks
    unsigned int nTracks = 100;

    // Set up RNG
    std::random_device rd;
    std::mt19937       gen(rd());

    // Set up constant B-Field
    ConstantBField bField(Vector3D(0., 0., 1.) * units::_T);

    // Set up Eigenstepper
    EigenStepper<ConstantBField> stepper(bField);

    // Set up propagator with void navigator
    Propagator<EigenStepper<ConstantBField>> propagator(stepper);

    // Create perigee surface
    std::shared_ptr<PerigeeSurface> perigeeSurface
        = Surface::makeShared<PerigeeSurface>(Vector3D(0., 0., 0.));

    // Create position of vertex and perigee surface
    double x = vXYDist(gen);
    double y = vXYDist(gen);
    double z = vZDist(gen);

    // Calculate d0 and z0 corresponding to vertex position
    double d0v = sqrt(x * x + y * y);
    double z0v = z;

    // Start constructing nTracks tracks in the following
    std::vector<BoundParameters> tracks;

    // Construct random track emerging from vicinity of vertex position
    // Vector to store track objects used for vertex fit
    for (unsigned int iTrack = 0; iTrack < nTracks; iTrack++) {
      // Construct positive or negative charge randomly
      double q = qDist(gen) < 0 ? -1. : 1.;

      // Construct random track parameters
      TrackParametersBase::ParVector_t paramVec;
      paramVec << d0v + d0Dist(gen), z0v + z0Dist(gen), phiDist(gen),
          thetaDist(gen), q / pTDist(gen);

      // Fill vector of track objects with simple covariance matrix
      std::unique_ptr<ActsSymMatrixD<5>> covMat
          = std::make_unique<ActsSymMatrixD<5>>();

      // Resolutions
      double resD0 = resIPDist(gen);
      double resZ0 = resIPDist(gen);
      double resPh = resAngDist(gen);
      double resTh = resAngDist(gen);
      double resQp = resQoPDist(gen);

      (*covMat) << resD0 * resD0, 0., 0., 0., 0., 0., resZ0 * resZ0, 0., 0., 0.,
          0., 0., resPh * resPh, 0., 0., 0., 0., 0., resTh * resTh, 0., 0., 0.,
          0., 0., resQp * resQp;
      tracks.push_back(
          BoundParameters(std::move(covMat), paramVec, perigeeSurface));
    }

    LinearizedTrackFactory<ConstantBField,
                           Propagator<EigenStepper<ConstantBField>>>::Config
        ltConfig(bField);
    LinearizedTrackFactory<ConstantBField,
                           Propagator<EigenStepper<ConstantBField>>>
        linFactory(ltConfig);

    ActsVectorD<5>    vec5Zero = ActsVectorD<5>::Zero();
    ActsSymMatrixD<5> mat5Zero = ActsSymMatrixD<5>::Zero();
    Vector3D          vec3Zero = Vector3D::Zero();
    ActsMatrixD<5, 3> mat53Zero = ActsMatrixD<5, 3>::Zero();

    for (const BoundParameters& parameters : tracks) {
      LinearizedTrack linTrack = linFactory.linearizeTrack(
          &parameters, Vector3D(0., 0., 0.), propagator);

      BOOST_CHECK_NE(linTrack.parametersAtPCA, vec5Zero);
      BOOST_CHECK_NE(linTrack.covarianceAtPCA, mat5Zero);
      BOOST_CHECK_EQUAL(linTrack.linearizationPoint, vec3Zero);
      BOOST_CHECK_NE(linTrack.positionJacobian, mat53Zero);
      BOOST_CHECK_NE(linTrack.momentumJacobian, mat53Zero);
      BOOST_CHECK_NE(linTrack.constantTerm, vec5Zero);
    }
  }

  BOOST_AUTO_TEST_CASE(linearized_track_factory_empty_test)
  {
    ConstantBField bField(Vector3D(0., 0., 1.) * units::_T);

    // Set up Eigenstepper
    EigenStepper<ConstantBField> stepper(bField);

    // Set up propagator with void navigator
    Propagator<EigenStepper<ConstantBField>> propagator(stepper);

    LinearizedTrackFactory<ConstantBField,
                           Propagator<EigenStepper<ConstantBField>>>::Config
        ltConfig(bField);
    LinearizedTrackFactory<ConstantBField,
                           Propagator<EigenStepper<ConstantBField>>>
        linFactory(ltConfig);

    ActsVectorD<5>    vec5Zero = ActsVectorD<5>::Zero();
    ActsSymMatrixD<5> mat5Zero = ActsSymMatrixD<5>::Zero();
    Vector3D          vec3Zero = Vector3D::Zero();
    ActsMatrixD<5, 3> mat53Zero = ActsMatrixD<5, 3>::Zero();

    LinearizedTrack linTrack
        = linFactory.linearizeTrack(nullptr, Vector3D(1., 2., 3.), propagator);

    BOOST_CHECK_EQUAL(linTrack.parametersAtPCA, vec5Zero);
    BOOST_CHECK_EQUAL(linTrack.covarianceAtPCA, mat5Zero);
    BOOST_CHECK_EQUAL(linTrack.linearizationPoint, vec3Zero);
    BOOST_CHECK_EQUAL(linTrack.positionJacobian, mat53Zero);
    BOOST_CHECK_EQUAL(linTrack.momentumJacobian, mat53Zero);
    BOOST_CHECK_EQUAL(linTrack.constantTerm, vec5Zero);
  }

}  // namespace Test
}  // namespace Acts

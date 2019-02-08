// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 Acts project team
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

namespace bdata = boost::unit_test::data;

namespace Acts {
namespace Test {

  struct InputTrack
  {
    InputTrack(const Acts::BoundParameters& params) : m_parameters(params) {}

    const Acts::BoundParameters&
    parameters() const
    {
      return m_parameters;
    }

  private:
    Acts::BoundParameters m_parameters;
  };

  const int ntests = 1;

  ///
  /// @brief Unit test for FullBilloirVertexFitter
  ///
  BOOST_DATA_TEST_CASE(
      linearized_track_factory_test,
      bdata::random((bdata::seed = 0,
                     bdata::distribution
                     = std::uniform_real_distribution<>(-0.1 * units::_mm,
                                                        0.1 * units::_mm)))
          ^ bdata::random(
                (bdata::seed = 1,
                 bdata::distribution
                 = std::uniform_real_distribution<>(-0.1 * units::_mm,
                                                    0.1 * units::_mm)))
          ^ bdata::random((bdata::seed = 2,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-20 * units::_mm,
                                                              20 * units::_mm)))
          ^ bdata::random((bdata::seed = 3,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-0.01, 0.01)))
          ^ bdata::random((bdata::seed = 4,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-0.01, 0.01)))
          ^ bdata::random((bdata::seed = 5,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-0.01, 0.01)))
          ^ bdata::random((bdata::seed = 6,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-0.01, 0.01)))
          ^ bdata::random((bdata::seed = 7,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-0.2, 0.2)))
          ^ bdata::random((bdata::seed = 8,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-0.2, 0.2)))
          ^ bdata::random((bdata::seed = 9,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-0.2, 0.2)))
          ^ bdata::random((bdata::seed = 10,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-0.2, 0.2)))
          ^ bdata::random(
                (bdata::seed = 11,
                 bdata::distribution
                 = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                    10. * units::_GeV)))
          ^ bdata::random(
                (bdata::seed = 12,
                 bdata::distribution
                 = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                    10. * units::_GeV)))
          ^ bdata::random(
                (bdata::seed = 13,
                 bdata::distribution
                 = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                    10. * units::_GeV)))
          ^ bdata::random(
                (bdata::seed = 14,
                 bdata::distribution
                 = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                    10. * units::_GeV)))
          ^ bdata::random((bdata::seed = 15,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::random((bdata::seed = 16,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::random((bdata::seed = 17,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::random((bdata::seed = 18,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::random((bdata::seed = 19,
                           bdata::distribution
                           = std::uniform_real_distribution<>(1.0, M_PI - 1.0)))
          ^ bdata::random((bdata::seed = 20,
                           bdata::distribution
                           = std::uniform_real_distribution<>(1.0, M_PI - 1.0)))
          ^ bdata::random((bdata::seed = 21,
                           bdata::distribution
                           = std::uniform_real_distribution<>(1.0, M_PI - 1.0)))
          ^ bdata::random((bdata::seed = 22,
                           bdata::distribution
                           = std::uniform_real_distribution<>(1.0, M_PI - 1.0)))
          ^ bdata::random((bdata::seed = 23,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-1, 1)))
          ^ bdata::random((bdata::seed = 24,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-1, 1)))
          ^ bdata::random((bdata::seed = 25,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-1, 1)))
          ^ bdata::random((bdata::seed = 26,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-1, 1)))
          ^ bdata::random(
                (bdata::seed = 27,
                 bdata::distribution
                 = std::uniform_real_distribution<>(0., 100. * units::_um)))
          ^ bdata::random(
                (bdata::seed = 28,
                 bdata::distribution
                 = std::uniform_real_distribution<>(0., 100. * units::_um)))
          ^ bdata::random((bdata::seed = 29,
                           bdata::distribution
                           = std::uniform_real_distribution<>(0., 0.1)))
          ^ bdata::random((bdata::seed = 30,
                           bdata::distribution
                           = std::uniform_real_distribution<>(0., 0.1)))
          ^ bdata::random((bdata::seed = 31,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-0.1, 0.1)))

          ^ bdata::xrange(ntests),
      x,
      y,
      z,
      d0_1,
      d0_2,
      d0_3,
      d0_4,
      z0_1,
      z0_2,
      z0_3,
      z0_4,
      pT_1,
      pT_2,
      pT_3,
      pT_4,
      phi_1,
      phi_2,
      phi_3,
      phi_4,
      theta_1,
      theta_2,
      theta_3,
      theta_4,
      q_1,
      q_2,
      q_3,
      q_4,
      res_d0,
      res_z0,
      res_ph,
      res_th,
      res_qp,

      index)
  {

    (void)index;

    // Store parameters for 4 tracks
    std::vector<double> d0Vec    = {d0_1, d0_2, d0_3, d0_4};
    std::vector<double> z0Vec    = {z0_1, z0_2, z0_3, z0_4};
    std::vector<double> pTVec    = {pT_1, pT_2, pT_3, pT_4};
    std::vector<double> phiVec   = {phi_1, phi_2, phi_3, phi_4};
    std::vector<double> thetaVec = {theta_1, theta_2, theta_3, theta_4};

    // Construct random vector of positive or negaitve charges
    std::vector<double> qTemp = {q_1, q_2, q_3, q_4};
    std::vector<int>    qVec;
    for (auto q : qTemp) {
      qVec.push_back(q < 0 ? -1. : 1.);
    }

    // Set up constant B-Field
    ConstantBField bField(Vector3D(0., 0., 1.) * units::_T);

    // Vector to store track objects used for vertex fit
    std::vector<InputTrack> tracks;

    // Create perigee surface
    std::shared_ptr<PerigeeSurface> perigeeSurface
        = Surface::makeShared<PerigeeSurface>(Vector3D(0., 0., 0.));

    // Calculate d0 and z0 corresponding to vertex position
    double d0_v = sqrt(x * x + y * y);
    double z0_v = z;

    // Start constructing 4 tracks in the following
    // Construct random track emerging from vicinity of vertex position
    for (unsigned int iTrack = 0; iTrack < d0Vec.size(); iTrack++) {
      TrackParametersBase::ParVector_t paramVec;
      paramVec << d0_v + d0Vec[iTrack], z0_v + z0Vec[iTrack], phiVec[iTrack],
          thetaVec[iTrack], ((double)qVec[iTrack]) / pTVec[iTrack];

      /// Fill vector of track objects with simple covariance matrix
      std::unique_ptr<ActsSymMatrixD<5>> covMat
          = std::make_unique<ActsSymMatrixD<5>>();
      (*covMat) << res_d0 * res_d0, 0., 0., 0., 0., 0., res_z0 * res_z0, 0., 0.,
          0., 0., 0., res_ph * res_ph, 0., 0., 0., 0., 0., res_th * res_th, 0.,
          0., 0., 0., 0., res_qp * res_qp;
      tracks.push_back(InputTrack(
          BoundParameters(std::move(covMat), paramVec, perigeeSurface)));
    }

    LinearizedTrackFactory<ConstantBField>::Config lt_config(bField);
    LinearizedTrackFactory<ConstantBField>         linFactory(lt_config);

    ActsVectorD<5>    vec5Zero = ActsVectorD<5>::Zero();
    ActsSymMatrixD<5> mat5Zero = ActsSymMatrixD<5>::Zero();
    Vector3D          vec3Zero = Vector3D::Zero();
    ActsMatrixD<5, 3> mat53Zero = ActsMatrixD<5, 3>::Zero();

    for (const InputTrack& track : tracks) {
      const auto&     trackParams = track.parameters();
      LinearizedTrack linTrack
          = linFactory.linearizeTrack(&trackParams, Vector3D(0., 0., 0.));

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
    LinearizedTrackFactory<ConstantBField>::Config lt_config(bField);
    LinearizedTrackFactory<ConstantBField>         linFactory(lt_config);

    ActsVectorD<5>    vec5Zero = ActsVectorD<5>::Zero();
    ActsSymMatrixD<5> mat5Zero = ActsSymMatrixD<5>::Zero();
    Vector3D          vec3Zero = Vector3D::Zero();
    ActsMatrixD<5, 3> mat53Zero = ActsMatrixD<5, 3>::Zero();

    LinearizedTrack linTrack
        = linFactory.linearizeTrack(nullptr, Vector3D(1., 2., 3.));

    BOOST_CHECK_EQUAL(linTrack.parametersAtPCA, vec5Zero);
    BOOST_CHECK_EQUAL(linTrack.covarianceAtPCA, mat5Zero);
    BOOST_CHECK_EQUAL(linTrack.linearizationPoint, vec3Zero);
    BOOST_CHECK_EQUAL(linTrack.positionJacobian, mat53Zero);
    BOOST_CHECK_EQUAL(linTrack.momentumJacobian, mat53Zero);
    BOOST_CHECK_EQUAL(linTrack.constantTerm, vec5Zero);
  }

}  // namespace Test
}  // namespace Acts

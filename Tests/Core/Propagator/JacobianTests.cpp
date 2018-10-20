// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE JacobianTests Tests
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/AtlasStepper.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  using BFieldType       = ConstantBField;
  using EigenStepperType = EigenStepper<BFieldType>;
  using AtlasStepperType = AtlasStepper<BFieldType>;

  /// Helper method to create a transform for a plane
  /// to mimic detector situations, the plane is roughly
  /// perpendicular to the track
  ///
  /// @param nnomal The nominal normal direction
  /// @param angleT Rotation around the norminal normal
  /// @param angleU Roation around the original U axis
  std::shared_ptr<Transform3D>
  createCylindricTransform(const Vector3D& nposition,
                           double          angleX,
                           double          angleY)
  {
    Transform3D ctransform;
    ctransform.setIdentity();
    ctransform.pretranslate(nposition);
    ctransform.prerotate(AngleAxis3D(angleX, Vector3D::UnitX()));
    ctransform.prerotate(AngleAxis3D(angleY, Vector3D::UnitY()));
    return std::make_shared<Transform3D>(ctransform);
  }

  /// Helper method to create a transform for a plane
  /// to mimic detector situations, the plane is roughly
  /// perpendicular to the track
  ///
  /// @param nnomal The nominal normal direction
  /// @param angleT Rotation around the norminal normal
  /// @param angleU Roation around the original U axis
  std::shared_ptr<Transform3D>
  createPlanarTransform(const Vector3D& nposition,
                        const Vector3D& nnormal,
                        double          angleT,
                        double          angleU)
  {
    // the rotation of the destination surface
    Vector3D T = nnormal.normalized();
    Vector3D U = std::abs(T.dot(Vector3D::UnitZ())) < 0.99
        ? Vector3D::UnitZ().cross(T).normalized()
        : Vector3D::UnitX().cross(T).normalized();
    Vector3D V = T.cross(U);
    // that's the plane curvilinear Rotation
    RotationMatrix3D curvilinearRotation;
    curvilinearRotation.col(0) = U;
    curvilinearRotation.col(1) = V;
    curvilinearRotation.col(2) = T;
    // curvilinear surfaces are boundless
    Transform3D ctransform{curvilinearRotation};
    ctransform.pretranslate(nposition);
    ctransform.prerotate(AngleAxis3D(angleT, T));
    ctransform.prerotate(AngleAxis3D(angleU, U));
    //
    return std::make_shared<Transform3D>(ctransform);
  }

  /// Helper method : convert into Acts matrix
  /// It takes the double array from AtlasStepper
  /// and transforms it into an ActsMatrixD
  ///
  /// @param P is the pointer to the array
  ///
  /// Translation is (for lookup)
  ///                   /dL0    /dL1    /dPhi   /dThe   /dCM
  /// X  ->P[0]  dX /   P[ 7]   P[14]   P[21]   P[28]   P[35]
  /// Y  ->P[1]  dY /   P[ 8]   P[15]   P[22]   P[29]   P[36]
  /// Z  ->P[2]  dZ /   P[ 9]   P[16]   P[23]   P[30]   P[37]
  /// Ax ->P[3]  dAx/   P[10]   P[17]   P[24]   P[31]   P[38]
  /// Ay ->P[4]  dAy/   P[11]   P[18]   P[25]   P[32]   P[39]
  /// Az ->P[5]  dAz/   P[12]   P[19]   P[26]   P[33]   P[40]
  /// CM ->P[6]  dCM/   P[13]   P[20]   P[27]   P[34]   P[41]

  ActsMatrixD<7, 5>
  convertToMatrix(const double* P)
  {
    // initialize to zero
    ActsMatrixD<7, 5> jMatrix = ActsMatrixD<7, 5>::Zero();
    for (size_t j = 0; j < 5; ++j) {
      for (size_t i = 0; i < 7; ++i) {
        size_t ijc = 7 + j * 7 + i;
        jMatrix(i, j) = P[ijc];
      }
    }
    return jMatrix;
  }

  /// Helper method : tests the jacobian to Global
  /// for a templated Parameters object
  ///
  /// @tparam Parameters the parameter type
  /// @param pars the parameter object
  template <typename Parameters>
  void
  testJacobianToGlobal(const Parameters& pars)
  {
    // Jacobian creation for Propagator/Steppers
    // a) ATLAS stepper
    AtlasStepperType::State astepState(pars);
    // b) Eigen stepper
    EigenStepperType::State estepState(pars);

    // create the matrices
    auto asMatrix = convertToMatrix(astepState.pVector);

    // cross comparison checks
    CHECK_CLOSE_OR_SMALL(asMatrix, estepState.jacToGlobal, 1e-6, 1e-9);
  }

  /// This tests the jacobian of local curvilinear -> global
  BOOST_AUTO_TEST_CASE(JacobianCurvilinearToGlobalTest)
  {
    ActsSymMatrixD<NGlobalPars> cov;
    cov << 10 * units::_mm, 0, 0, 0, 0, 0, 10 * units::_mm, 0, 0, 0, 0, 0, 0.1,
        0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 1. / (10 * units::_GeV);
    auto covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);

    // Let's create a surface somewhere in space
    Vector3D     pos(341., 412., 93.);
    Vector3D     mom(1.2, 8.3, 0.45);
    const double q = 1;

    // Create curvilinear parameters
    CurvilinearParameters curvilinear(std::move(covPtr), pos, mom, q);

    // run the test
    testJacobianToGlobal(curvilinear);
  }

  /// This tests the jacobian of local cylinder -> global
  BOOST_AUTO_TEST_CASE(JacobianCylinderToGlobalTest)
  {
    // the cylinder transform and surface
    auto cTransform = createCylindricTransform({10., -5., 0.}, 0.004, 0.03);
    auto cSurface
        = Surface::makeShared<CylinderSurface>(cTransform, 200., 1000.);

    ActsSymMatrixD<NGlobalPars> cov;
    cov << 10 * units::_mm, 0, 0, 0, 0, 0, 10 * units::_mm, 0, 0, 0, 0, 0, 0.1,
        0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 1. / (10 * units::_GeV);
    auto covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);

    ActsVectorD<NGlobalPars> pars;
    pars << 182.34, -82., 0.134, 0.85, 1. / (100 * units::_GeV);

    BoundParameters atCylinder(std::move(covPtr), std::move(pars), cSurface);

    // run the test
    testJacobianToGlobal(atCylinder);
  }

  /// This tests the jacobian of local disc -> global
  BOOST_AUTO_TEST_CASE(JacobianDiscToGlobalTest)
  {

    // the disc transform and surface
    auto dTransform = createPlanarTransform(
        {10., -5., 0.}, Vector3D(0.23, 0.07, 1.).normalized(), 0.004, 0.03);
    auto dSurface = Surface::makeShared<DiscSurface>(dTransform, 200., 1000.);

    ActsSymMatrixD<NGlobalPars> cov;
    cov << 10 * units::_mm, 0, 0, 0, 0, 0, 10 * units::_mm, 0, 0, 0, 0, 0, 0.1,
        0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 1. / (10 * units::_GeV);
    auto covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);

    ActsVectorD<NGlobalPars> pars;
    pars << 192.34, 1.823, 0.734, 0.235, 1. / (100 * units::_GeV);

    BoundParameters atDisc(std::move(covPtr), std::move(pars), dSurface);

    // run the test
    testJacobianToGlobal(atDisc);
  }

  /// This tests the jacobian of local plane -> global
  BOOST_AUTO_TEST_CASE(JacobianPlaneToGlobalTest)
  {
    // Let's create a surface somewhere in space
    Vector3D sPosition(3421., 112., 893.);
    Vector3D sNormal = Vector3D(1.2, -0.3, 0.05).normalized();

    // Create a surface & parameters with covariance on the surface
    auto pSurface = Surface::makeShared<PlaneSurface>(sPosition, sNormal);

    ActsSymMatrixD<NGlobalPars> cov;
    cov << 10 * units::_mm, 0, 0, 0, 0, 0, 10 * units::_mm, 0, 0, 0, 0, 0, 0.1,
        0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 1. / (10 * units::_GeV);
    auto covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);

    ActsVectorD<NGlobalPars> pars;
    pars << 12.34, -8722., 2.134, 0.85, 1. / (100 * units::_GeV);

    BoundParameters atPlane(std::move(covPtr), std::move(pars), pSurface);

    // run the test
    testJacobianToGlobal(atPlane);
  }

  /// This tests the jacobian of local perigee -> global
  BOOST_AUTO_TEST_CASE(JacobianPerigeeToGlobalTest)
  {

    // Create a surface & parameters with covariance on the surface
    auto pSurface = Surface::makeShared<PerigeeSurface>(Vector3D({0., 0., 0.}));

    ActsSymMatrixD<NGlobalPars> cov;
    cov << 10 * units::_mm, 0, 0, 0, 0, 0, 10 * units::_mm, 0, 0, 0, 0, 0, 0.1,
        0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 1. / (10 * units::_GeV);
    auto covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);

    ActsVectorD<NGlobalPars> pars;
    pars << -3.34, -822., -0.734, 0.85, 1. / (100 * units::_GeV);

    BoundParameters perigee(std::move(covPtr), std::move(pars), pSurface);

    // run the test
    testJacobianToGlobal(perigee);
  }

  /// This tests the jacobian of local straw -> global
  BOOST_AUTO_TEST_CASE(JacobianStrawToGlobalTest)
  {
    // Create a surface & parameters with covariance on the surface
    auto sTransform = createCylindricTransform({1019., -52., 382.}, 0.4, -0.3);
    auto sSurface   = Surface::makeShared<StrawSurface>(sTransform, 10., 1000.);

    ActsSymMatrixD<NGlobalPars> cov;
    cov << 10 * units::_mm, 0, 0, 0, 0, 0, 10 * units::_mm, 0, 0, 0, 0, 0, 0.1,
        0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 1. / (10 * units::_GeV);
    auto covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);

    ActsVectorD<NGlobalPars> pars;
    pars << -8.34, 812., 0.734, 0.25, 1. / (100 * units::_GeV);

    BoundParameters atStraw(std::move(covPtr), std::move(pars), sSurface);

    // run the test
    testJacobianToGlobal(atStraw);
  }

}  // namespace Test
}  // namespace Acts

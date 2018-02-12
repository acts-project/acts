
// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE JacobianTests Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/detail/RungeKuttaUtils.hpp"
#include "ACTS/MagneticField/ConstantBField.hpp"
#include "ACTS/Propagator/AtlasStepper.hpp"
#include "ACTS/Propagator/EigenStepper.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Surfaces/PerigeeSurface.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  typedef ConstantBField            BField_type;
  typedef EigenStepper<BField_type> EigenStepper_type;
  typedef AtlasStepper<BField_type> AtlasStepper_type;

  /// Helper method to create a transform for a plane
  /// to mimic detector situations, the plane is roughtly
  /// perpenticular to the track
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
  /// to mimic detector situations, the plane is roughtly
  /// perpenticular to the track
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
    // that's the plane curinilear Rotation
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
  /// It takes the double array from AtlasStepper and RungeKuttaUtils
  /// and transforms it into an ActsMatrixD
  ///
  /// @param P is the pointer to the array
  ActsMatrixD<7, 5>
  convertToMatrix(const double* P)
  {
    // initialize to zero
    ActsMatrixD<7, 5> jMatrix = ActsMatrixD<7, 5>::Zero();
    for (size_t j = 0; j < 5; ++j)
      for (size_t i = 0; i < 7; ++i) {
        size_t ijc = 7 + j * 7 + i;
        jMatrix(i, j) = P[ijc];
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
    //
    // a) Original ATLAS code - using RungeKuttaUtils
    const size_t PDIM = 64;
    double       P[PDIM];
    // initialize
    RungeKuttaUtils rkUtils;
    auto            tfl = rkUtils.transformLocalToGlobal(true, pars, P);
    BOOST_CHECK(tfl);
    // b) ATLAS stepper
    AtlasStepper_type::Cache asCache(pars);
    // c) Eigen stepper
    EigenStepper_type::Cache esCache(pars);

    // create the matrices
    auto rkMatrix = convertToMatrix(P);
    auto asMatrix = convertToMatrix(asCache.pVector);

    // cross comparison checks
    BOOST_CHECK(rkMatrix.isApprox(asMatrix));
    BOOST_CHECK(rkMatrix.isApprox(esCache.jacobian));
    BOOST_CHECK(asMatrix.isApprox(esCache.jacobian));
  }

  /// This tests the jacobian of local curvilinear -> global
  BOOST_AUTO_TEST_CASE(JacobianCurvilinearToGlobalTest)
  {

    ActsSymMatrixD<NGlobalPars> cov;
    cov << 10 * units::_mm, 0, 0, 0, 0, 0, 10 * units::_mm, 0, 0, 0, 0, 0, 0.1,
        0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 1. / (10 * units::_GeV);
    auto cov_ptr = std::make_unique<const ActsSymMatrixD<5>>(cov);

    // Let's create a surface somewhere in space
    Vector3D     pos(341., 412., 93.);
    Vector3D     mom(1.2, 8.3, 0.45);
    const double q = 1;

    // Create curilinear parameters
    CurvilinearParameters curvilinear(std::move(cov_ptr), pos, mom, q);

    // run the test
    testJacobianToGlobal(curvilinear);
  }

  /// This tests the jacobian of local cylinder -> global
  BOOST_AUTO_TEST_CASE(JacobianCylinderToGlobalTest)
  {

    // the cylinder transform and surface
    auto cTransform = createCylindricTransform({10., -5., 0.}, 0.004, 0.03);
    CylinderSurface cSurface(cTransform, 200., 1000.);

    ActsSymMatrixD<NGlobalPars> cov;
    cov << 10 * units::_mm, 0, 0, 0, 0, 0, 10 * units::_mm, 0, 0, 0, 0, 0, 0.1,
        0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 1. / (10 * units::_GeV);
    auto cov_ptr = std::make_unique<const ActsSymMatrixD<5>>(cov);

    ActsVectorD<NGlobalPars> pars;
    pars << 182.34, -82., 0.134, 0.85, 1. / (100 * units::_GeV);

    BoundParameters atCylinder(std::move(cov_ptr), std::move(pars), cSurface);

    // run the test
    testJacobianToGlobal(atCylinder);
  }

  /// This tests the jacobian of local plane -> global
  BOOST_AUTO_TEST_CASE(JacobianPlaneToGlobalTest)
  {
    // Let's create a surface somewhere in space
    Vector3D sPosition(3421., 112., 893.);
    Vector3D sNormal = Vector3D(1.2, -0.3, 0.05).normalized();

    // Create a surface & parameters with covariance on the surface
    PlaneSurface pSurface(sPosition, sNormal);

    ActsSymMatrixD<NGlobalPars> cov;
    cov << 10 * units::_mm, 0, 0, 0, 0, 0, 10 * units::_mm, 0, 0, 0, 0, 0, 0.1,
        0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 1. / (10 * units::_GeV);
    auto cov_ptr = std::make_unique<const ActsSymMatrixD<5>>(cov);

    ActsVectorD<NGlobalPars> pars;
    pars << 12.34, -8722., 2.134, 0.85, 1. / (100 * units::_GeV);

    BoundParameters atPlane(std::move(cov_ptr), std::move(pars), pSurface);

    // run the test
    testJacobianToGlobal(atPlane);
  }

  /// This tests the jacobian of local line -> global
  BOOST_AUTO_TEST_CASE(JacobianPerigeeToGlobalTest)
  {

    // Create a surface & parameters with covariance on the surface
    PerigeeSurface pSurface({0., 0., 0.});

    ActsSymMatrixD<NGlobalPars> cov;
    cov << 10 * units::_mm, 0, 0, 0, 0, 0, 10 * units::_mm, 0, 0, 0, 0, 0, 0.1,
        0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 1. / (10 * units::_GeV);
    auto cov_ptr = std::make_unique<const ActsSymMatrixD<5>>(cov);

    ActsVectorD<NGlobalPars> pars;
    pars << -3.34, -822., -0.734, 0.85, 1. / (100 * units::_GeV);

    BoundParameters perigee(std::move(cov_ptr), std::move(pars), pSurface);

    // run the test
    testJacobianToGlobal(perigee);
  }

}  // namespace Test
}  // namespace Acts

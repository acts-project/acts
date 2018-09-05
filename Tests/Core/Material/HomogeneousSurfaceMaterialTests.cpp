// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE HomogeneousSurfaceMaterial Tests
#include <boost/test/included/unit_test.hpp>
#include <climits>
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"

namespace Acts {

namespace Test {

  /// Test the constructors
  BOOST_AUTO_TEST_CASE(HomogeneousSurfaceMaterial_construction_test)
  {
    // construct the material properties from argumnets
    MaterialProperties mp(1., 2., 3., 4., 5., 0.1);

    // Constructor from arguments
    HomogeneousSurfaceMaterial hsm(mp, 1.);
    // Copy constructor
    HomogeneousSurfaceMaterial hsmCopy(hsm);
    // Test equality of the copy
    BOOST_TEST(hsm == hsmCopy);
    // Copy move constructor
    HomogeneousSurfaceMaterial hsmCopyMoved(std::move(hsmCopy));
    // Test equality of the copy
    BOOST_TEST(hsm == hsmCopyMoved);
    // Assignment constructor
    HomogeneousSurfaceMaterial hsmAssigned = hsm;
    // Test equality of the asignment
    BOOST_TEST(hsm == hsmAssigned);
    // Assignment move constructor
    HomogeneousSurfaceMaterial hsmAssignedMoved(std::move(hsmAssigned));
    // Test equality of the copy
    BOOST_TEST(hsm == hsmAssignedMoved);
  }

  // Test the Scaling
  BOOST_AUTO_TEST_CASE(HomogeneousSurfaceMaterial_scaling_test)
  {

    // Construct the material properties from argumnets
    MaterialProperties mat(1., 2., 3., 4., 5., 0.1);
    MaterialProperties matHalf = mat;
    matHalf *= 0.5;

    HomogeneousSurfaceMaterial hsm(mat, 1.);
    hsm *= 0.5;

    auto matBin = hsm.materialProperties(0, 0);

    BOOST_TEST(matBin == matHalf);
    BOOST_TEST(matBin != mat);
  }

  // Test the Access
  BOOST_AUTO_TEST_CASE(HomogeneousSurfaceMaterial_access_test)
  {
    // construct the material properties from argumnets
    MaterialProperties mat(1., 2., 3., 4., 5., 0.1);
    MaterialProperties matHalf = mat;
    matHalf *= 0.5;

    MaterialProperties vacuum = MaterialProperties();

    // Constructor from arguments
    HomogeneousSurfaceMaterial hsmfwd(mat, 1.);
    HomogeneousSurfaceMaterial hsmhalf(mat, 0.5);
    HomogeneousSurfaceMaterial hsmbwd(mat, 0.);

    auto mat2d  = hsmfwd.materialProperties(Vector2D{0., 0.});
    auto mat3d  = hsmfwd.materialProperties(Vector3D{0., 0., 0.});
    auto matbin = hsmfwd.materialProperties(0, 0);

    // Test equality of the copy
    BOOST_TEST(mat == mat2d);
    BOOST_TEST(mat == mat3d);
    BOOST_TEST(mat == matbin);

    NavigationDirection fDir = forward;
    NavigationDirection bDir = backward;

    MaterialUpdateStage pre  = preUpdate;
    MaterialUpdateStage full = fullUpdate;
    MaterialUpdateStage post = postUpdate;

    // (a) Forward factor material test
    BOOST_TEST(hsmfwd.factor(fDir, full) == 1.);
    BOOST_TEST(hsmfwd.factor(fDir, pre) == 0.);
    BOOST_TEST(hsmfwd.factor(fDir, post) == 1.);

    BOOST_TEST(hsmfwd.factor(bDir, full) == 1.);
    BOOST_TEST(hsmfwd.factor(bDir, pre) == 1.);
    BOOST_TEST(hsmfwd.factor(bDir, post) == 0.);

    auto matFwdFull
        = hsmfwd.materialProperties(Vector3D{0., 0., 0.}, fDir, full);
    auto matBwdFull
        = hsmfwd.materialProperties(Vector3D{0., 0., 0.}, bDir, full);

    auto matFwdPost
        = hsmfwd.materialProperties(Vector3D{0., 0., 0.}, fDir, post);
    auto matBwdPost
        = hsmfwd.materialProperties(Vector3D{0., 0., 0.}, bDir, post);

    auto matFwdPre = hsmfwd.materialProperties(Vector3D{0., 0., 0.}, fDir, pre);
    auto matBwdPre = hsmfwd.materialProperties(Vector3D{0., 0., 0.}, bDir, pre);

    BOOST_TEST(mat == matFwdFull);
    BOOST_TEST(mat == matBwdFull);

    BOOST_TEST(mat == matFwdPost);
    BOOST_TEST(vacuum == matBwdPost);

    BOOST_TEST(vacuum == matFwdPre);
    BOOST_TEST(mat == matBwdPre);

    // (b) Split factor material test
    BOOST_TEST(hsmhalf.factor(fDir, full) == 1.);
    BOOST_TEST(hsmhalf.factor(fDir, pre) == 0.5);
    BOOST_TEST(hsmhalf.factor(fDir, post) == 0.5);

    BOOST_TEST(hsmhalf.factor(bDir, full) == 1.);
    BOOST_TEST(hsmhalf.factor(bDir, pre) == 0.5);
    BOOST_TEST(hsmhalf.factor(bDir, post) == 0.5);

    matFwdFull = hsmhalf.materialProperties(Vector3D{0., 0., 0.}, fDir, full);
    matBwdFull = hsmhalf.materialProperties(Vector3D{0., 0., 0.}, bDir, full);

    matFwdPost = hsmhalf.materialProperties(Vector3D{0., 0., 0.}, fDir, post);
    matBwdPost = hsmhalf.materialProperties(Vector3D{0., 0., 0.}, bDir, post);

    matFwdPre = hsmhalf.materialProperties(Vector3D{0., 0., 0.}, fDir, pre);
    matBwdPre = hsmhalf.materialProperties(Vector3D{0., 0., 0.}, bDir, pre);

    BOOST_TEST(mat == matFwdFull);
    BOOST_TEST(mat == matBwdFull);

    BOOST_TEST(matHalf == matFwdPost);
    BOOST_TEST(matHalf == matBwdPost);

    BOOST_TEST(matHalf == matFwdPre);
    BOOST_TEST(matHalf == matBwdPre);

    // c) Forward factor material test
    BOOST_TEST(hsmbwd.factor(fDir, full) == 1.);
    BOOST_TEST(hsmbwd.factor(fDir, pre) == 1.);
    BOOST_TEST(hsmbwd.factor(fDir, post) == 0.);

    BOOST_TEST(hsmbwd.factor(bDir, full) == 1.);
    BOOST_TEST(hsmbwd.factor(bDir, pre) == 0.);
    BOOST_TEST(hsmbwd.factor(bDir, post) == 1.);

    matFwdFull = hsmbwd.materialProperties(Vector3D{0., 0., 0.}, fDir, full);
    matBwdFull = hsmbwd.materialProperties(Vector3D{0., 0., 0.}, bDir, full);

    matFwdPost = hsmbwd.materialProperties(Vector3D{0., 0., 0.}, fDir, post);
    matBwdPost = hsmbwd.materialProperties(Vector3D{0., 0., 0.}, bDir, post);

    matFwdPre = hsmbwd.materialProperties(Vector3D{0., 0., 0.}, fDir, pre);
    matBwdPre = hsmbwd.materialProperties(Vector3D{0., 0., 0.}, bDir, pre);

    BOOST_TEST(mat == matFwdFull);
    BOOST_TEST(mat == matBwdFull);

    BOOST_TEST(vacuum == matFwdPost);
    BOOST_TEST(mat == matBwdPost);

    BOOST_TEST(mat == matFwdPre);
    BOOST_TEST(vacuum == matBwdPre);
  }
}
}

// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file ConstantBField_tests.cpp

#define BOOST_TEST_MODULE Constant magnetic field tests
#include <boost/test/included/unit_test.hpp>

#include <boost/test/data/test_case.hpp>
#include "ACTS/MagneticField/ConstantBField.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

   
  /// @brief unit test for construction of constant magnetic field
  ///
  /// Tests the correct behavior and consistency of
  /// -# ConstantBField::ConstantBField(double Bx,double By,double Bz)
  /// -# ConstantBField::ConstantBField(Vector3D B)
  /// -# ConstantBField::getField(const double* xyz, double* B) const
  /// -# ConstantBField::getField(const Vector3D& pos) const
  BOOST_DATA_TEST_CASE(ConstantBField_components,
                       bdata::random(-2. * units::_T, 2. * units::_T)
                           ^ bdata::random(-1. * units::_T, 4. * units::_T)
                           ^ bdata::random(0. * units::_T, 10. * units::_T)
                           ^ bdata::random(-10. * units::_m, 10. * units::_m)
                           ^ bdata::random(-10. * units::_m, 10. * units::_m)
                           ^ bdata::random(-10. * units::_m, 10. * units::_m)
                           ^ bdata::xrange(10),
                       x,
                       y,
                       z,
                       bx,
                       by,
                       bz,
                       index)
  {
    BOOST_TEST_CONTEXT("Eigen interface")
    {
      const Vector3D       Btrue(bx, by, bz);
      const Vector3D       pos(x, y, z);
      const ConstantBField BField(Btrue);

      BOOST_TEST(Btrue == BField.getField(pos));
      BOOST_TEST(Btrue == BField.getField(Vector3D(0, 0, 0)));
      BOOST_TEST(Btrue == BField.getField(-2 * pos));
    }

    BOOST_TEST_CONTEXT("C-array interface")
    {
      const double         Btrue[3] = {bx, by, bz};
      const double         pos[3]   = {x, y, z};
      const ConstantBField BField(bx, by, bz);
      double               Btest[3] = {0, 0, 0};

      BField.getField(pos, Btest);
      BOOST_TEST(Btrue[0] == Btest[0]);
      BOOST_TEST(Btrue[1] == Btest[1]);
      BOOST_TEST(Btrue[2] == Btest[2]);

      BField.getField(0, Btest);
      BOOST_TEST(Btrue[0] == Btest[0]);
      BOOST_TEST(Btrue[1] == Btest[1]);
      BOOST_TEST(Btrue[2] == Btest[2]);
    }

    BOOST_TEST_CONTEXT("Eigen initialised - C-array retrieved")
    {
      const Vector3D       Btrue(bx, by, bz);
      const ConstantBField BField(Btrue);
      const double         pos[3]   = {x, y, z};
      double               Btest[3] = {0, 0, 0};

      BField.getField(pos, Btest);
      BOOST_TEST(Btrue(0) == Btest[0]);
      BOOST_TEST(Btrue(1) == Btest[1]);
      BOOST_TEST(Btrue(2) == Btest[2]);

      BField.getField(0, Btest);
      BOOST_TEST(Btrue(0) == Btest[0]);
      BOOST_TEST(Btrue(1) == Btest[1]);
      BOOST_TEST(Btrue(2) == Btest[2]);
    }

    BOOST_TEST_CONTEXT("C-array initialised - Eigen retrieved")
    {
      const ConstantBField BField(bx, by, bz);
      const Vector3D       Btrue(bx, by, bz);
      const Vector3D       pos(x, y, z);

      BOOST_TEST(Btrue == BField.getField(pos));
      BOOST_TEST(Btrue == BField.getField(Vector3D(0, 0, 0)));
      BOOST_TEST(Btrue == BField.getField(-2 * pos));
    }
  }

  /// @brief unit test for update of constant magnetic field
  ///
  /// Tests the correct behavior and consistency of
  /// -# ConstantBField::setField(double Bx, double By, double Bz)
  /// -# ConstantBField::setField(const Vector3D& B)
  /// -# ConstantBField::getField(const double* xyz, double* B) const
  /// -# ConstantBField::getField(const Vector3D& pos) const
  BOOST_DATA_TEST_CASE(ConstantBField_update,
                       bdata::random(-2. * units::_T, 2. * units::_T)
                           ^ bdata::random(-1. * units::_T, 4. * units::_T)
                           ^ bdata::random(0. * units::_T, 10. * units::_T)
                           ^ bdata::random(-10. * units::_m, 10. * units::_m)
                           ^ bdata::random(-10. * units::_m, 10. * units::_m)
                           ^ bdata::random(-10. * units::_m, 10. * units::_m)
                           ^ bdata::xrange(10),
                       x,
                       y,
                       z,
                       bx,
                       by,
                       bz,
                       index)
  {
    ConstantBField BField(0, 0, 0);

    BOOST_TEST_CONTEXT("Eigen interface")
    {
      const Vector3D Btrue(bx, by, bz);
      const Vector3D pos(x, y, z);
      BField.setField(Btrue);

      BOOST_TEST(Btrue == BField.getField(pos));
      BOOST_TEST(Btrue == BField.getField(Vector3D(0, 0, 0)));
      BOOST_TEST(Btrue == BField.getField(-2 * pos));
    }

    BOOST_TEST_CONTEXT("C-array interface")
    {
      const double pos[3]   = {x, y, z};
      double       Btest[3] = {0, 0, 0};
      BField.setField(bx, by, bz);

      BField.getField(pos, Btest);
      BOOST_TEST(bx == Btest[0]);
      BOOST_TEST(by == Btest[1]);
      BOOST_TEST(bz == Btest[2]);

      BField.getField(0, Btest);
      BOOST_TEST(bx == Btest[0]);
      BOOST_TEST(by == Btest[1]);
      BOOST_TEST(bz == Btest[2]);
    }

    BOOST_TEST_CONTEXT("Eigen initialised - C-array retrieved")
    {
      const Vector3D Btrue(bx, by, bz);
      const double   pos[3]   = {x, y, z};
      double         Btest[3] = {0, 0, 0};
      BField.setField(Btrue);

      BField.getField(pos, Btest);
      BOOST_TEST(Btrue(0) == Btest[0]);
      BOOST_TEST(Btrue(1) == Btest[1]);
      BOOST_TEST(Btrue(2) == Btest[2]);

      BField.getField(0, Btest);
      BOOST_TEST(Btrue(0) == Btest[0]);
      BOOST_TEST(Btrue(1) == Btest[1]);
      BOOST_TEST(Btrue(2) == Btest[2]);
    }

    BOOST_TEST_CONTEXT("C-array initialised - Eigen retrieved")
    {
      const Vector3D Btrue(bx, by, bz);
      const Vector3D pos(x, y, z);
      BField.setField(bx, by, bz);

      BOOST_TEST(Btrue == BField.getField(pos));
      BOOST_TEST(Btrue == BField.getField(Vector3D(0, 0, 0)));
      BOOST_TEST(Btrue == BField.getField(-2 * pos));
    }
  }
}  // end of namespace Test

}  // end of namespace Acts

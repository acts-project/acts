// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <boost/test/included/unit_test.hpp>
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
namespace Test {

  template <typename T, int rows, int cols>
  void checkCloseImpl(const ActsMatrix<T, rows, cols>& m1,
                      const ActsMatrix<T, rows, cols>& m2)
  {
    const int m_size = rows * cols;
    for (int i = 0; i < m_size; ++i) {
      BOOST_CHECK_CLOSE_FRACTION(m1(i), m2(i), 1e-6);
    }
  }

  void checkCloseVec3D(const Vector3D& v1, const Vector3D& v2) {
    checkCloseImpl(v1, v2);
  }

  void checkCloseRM3D(const RotationMatrix3D& m1, const RotationMatrix3D& m2) {
    checkCloseImpl(m1, m2);
  }

  template <typename Parameter>
  void
  consistencyCheck(const Parameter& pars,
                   const Vector3D&  position,
                   const Vector3D&  momentum,
                   double           charge,
                   std::array<double, 5> values)
  {
    // check parameter vector
    BOOST_CHECK_CLOSE(
        pars.parameters()[eLOC_0], values[0], s_onSurfaceTolerance);
    BOOST_CHECK_CLOSE(
        pars.parameters()[eLOC_1], values[1], s_onSurfaceTolerance);
    BOOST_CHECK_CLOSE(pars.parameters()[ePHI], values[2], s_onSurfaceTolerance);
    BOOST_CHECK_CLOSE(
        pars.parameters()[eTHETA], values[3], s_onSurfaceTolerance);
    BOOST_CHECK_CLOSE(pars.parameters()[eQOP], values[4], s_onSurfaceTolerance);
    // check global parameters
    BOOST_CHECK(pars.position().isApprox(position));
    BOOST_CHECK(pars.momentum().isApprox(momentum));
    BOOST_CHECK_EQUAL(pars.charge(), charge);
  }
}
}
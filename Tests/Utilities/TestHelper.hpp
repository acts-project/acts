// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
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

  void checkCloseVec2D(const Vector2D& v1, const Vector2D& v2) {
    checkCloseImpl(v1, v2);
  }

  void checkCloseVec3D(const Vector3D& v1, const Vector3D& v2) {
    checkCloseImpl(v1, v2);
  }

  void checkCloseRM3D(const RotationMatrix3D& m1, const RotationMatrix3D& m2) {
    checkCloseImpl(m1, m2);
  }

}
}
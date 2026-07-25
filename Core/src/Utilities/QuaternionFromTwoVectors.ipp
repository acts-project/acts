// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/detail/QuaternionFromTwoVectors.hpp"

namespace Acts::detail {

template <typename T>
Eigen::Quaternion<T> quaternionFromTwoVectors(const Eigen::Matrix<T, 3, 1>& a,
                                              const Eigen::Matrix<T, 3, 1>& b) {
  return Eigen::Quaternion<T>().setFromTwoVectors(a, b);
}

}  // namespace Acts::detail

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "QuaternionFromTwoVectors.ipp"

// GCC (>=12) emits a spurious -Wmaybe-uninitialized from inside Eigen's SVD
// code used by setFromTwoVectors; suppress it for this translation unit.
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ >= 12
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

namespace Acts::detail {

template Eigen::Quaternion<double> quaternionFromTwoVectors<double>(
    const Eigen::Matrix<double, 3, 1>&, const Eigen::Matrix<double, 3, 1>&);

}  // namespace Acts::detail

#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ >= 12
#pragma GCC diagnostic pop
#endif

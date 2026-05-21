// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <cstddef>
#include <span>

#include <Eigen/Core>

namespace Acts::detail {

template <typename T, std::size_t N>
std::array<T, N> stdArrayCopy(const std::span<const T, N> arr) {
  std::array<T, N> result{};
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = arr[i];
  }
  return result;
}

template <typename T, int N>
std::array<T, N> stdArrayCopy(const Eigen::Vector<T, N>& arr)
  requires(N > 0)
{
  std::array<T, N> result{};
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = arr[i];
  }
  return result;
}

template <typename T, std::size_t N>
Eigen::Vector<T, N> stdArrayToEigen(const std::array<T, N>& arr) {
  Eigen::Vector<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = arr[i];
  }
  return result;
}

template <typename T, std::size_t N>
std::array<T, N> stdArrayAddScaled(const std::array<T, N>& a,
                                   const std::array<T, N>& b, const T scale) {
  std::array<T, N> result{};
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = a[i] + b[i] * scale;
  }
  return result;
}

template <typename T, std::size_t N>
T stdArrayDot(const std::array<T, N>& a, const std::array<T, N>& b) {
  T result = 0;
  for (std::size_t i = 0; i < N; ++i) {
    result += a[i] * b[i];
  }
  return result;
}

template <typename T>
std::array<T, 3> stdArrayCross(const std::array<T, 3>& a,
                               const std::array<T, 3>& b) {
  std::array<T, 3> result{};
  result[0] = a[1] * b[2] - a[2] * b[1];
  result[1] = a[2] * b[0] - a[0] * b[2];
  result[2] = a[0] * b[1] - a[1] * b[0];
  return result;
}

}  // namespace Acts::detail

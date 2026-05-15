// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <span>

namespace Acts::detail {

template <typename T, std::size_t N>
std::span<T, N> cmathMap(std::array<T, N>& arr) {
  return std::span<T, N>(arr.data(), arr.size());
}

template <typename T, std::size_t N>
std::span<const T, N> cmathMap(const std::array<T, N>& arr) {
  return std::span<const T, N>(arr.data(), arr.size());
}

template <typename T, std::size_t N>
void cmathCopy(const std::span<const T, N> input,
               const std::span<T, N> result) {
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = input[i];
  }
}

template <typename T, std::size_t N>
std::array<T, N> cmathCopy(const std::span<const T, N> input) {
  std::array<T, N> result{};
  cmathCopy(input, cmathMap(result));
  return result;
}

template <typename T, std::size_t N>
void cmathAddScaled(const std::span<const T, N> a,
                    const std::span<const T, N> b, const T scale,
                    const std::span<T, N> result) {
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = a[i] + b[i] * scale;
  }
}

template <typename T, std::size_t N>
std::array<T, N> cmathAddScaled(const std::span<const T, N> a,
                                const std::span<const T, N> b, const T scale) {
  std::array<T, N> result{};
  cmathAddScaled(a, b, scale, cmathMap(result));
  return result;
}

template <typename T, std::size_t N>
T cmathDot(const std::span<const T, N> a, const std::span<const T, N> b) {
  T result = 0;
  for (std::size_t i = 0; i < N; ++i) {
    result += a[i] * b[i];
  }
  return result;
}

template <typename T>
void cmathCross(const std::span<const T, 3> a, const std::span<const T, 3> b,
                const std::span<T, 3> result) {
  result[0] = a[1] * b[2] - a[2] * b[1];
  result[1] = a[2] * b[0] - a[0] * b[2];
  result[2] = a[0] * b[1] - a[1] * b[0];
}

template <typename T>
std::array<T, 3> cmathCross(const std::span<const T, 3> a,
                            const std::span<const T, 3> b) {
  std::array<T, 3> result{};
  cmathCross(a, b, cmathMap(result));
  return result;
}

}  // namespace Acts::detail

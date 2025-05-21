// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <array>

namespace Acts {
///  @brief creates an array of type T and size N and assigns all elements to the parsed
///         default value def_val. This basically allows for inline construction
///         and initialization of class-member arrays
///  @param def_val: Default value to assign
template <typename T, std::size_t N>
constexpr std::array<T, N> filledArray(const T& def_val) {
  std::array<T, N> arr{};
  arr.fill(def_val);
  return arr;
}
}  // namespace Acts

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>

namespace Acts {

template <typename T>
inline auto square(T x) {
  return x * x;
}

template <typename... T>
inline auto hypot2(T... args) {
  return (square(args) + ...);
}

template <typename... T>
inline auto hypot(T... args) {
  return std::sqrt(hypot2(args...));
}

}  // namespace Acts

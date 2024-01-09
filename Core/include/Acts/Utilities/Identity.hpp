// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#if __cplusplus >= 202002L

#include <functional>

namespace Acts {
using Identity = std::identity;
}

#else

#include <utility>

namespace Acts {

/// @brief Function object which maps a value to itself by perfect forwarding
/// This is a backport of C++20's std::identity
struct Identity {
  template <typename T>
  constexpr auto operator()(T &&v) const {
    return std::forward<T>(v);
  }
};

}  // namespace Acts

#endif

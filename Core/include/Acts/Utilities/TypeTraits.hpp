// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <type_traits>

namespace Acts {
/// @brief Type trait that adds const qualifier to a type based on a boolean condition
/// @tparam C Boolean condition determining whether to add const
/// @tparam T Type to potentially make const
template <bool C, typename T>
using const_if_t = std::conditional_t<C, const T, T>;
}  // namespace Acts

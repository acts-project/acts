// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// This workround does not work on libc++. To detect libc++, we include
// one STL header and then check if _LIBCPP_VERSION is defined.

#include <any>
#include <type_traits>

// only if clang and libstdc++
#if !defined(_LIBCPP_VERSION) && defined(__clang__)

// Workaround for building on clang+libstdc++
// See https://gitlab.cern.ch/atlas/atlasexternals/merge_requests/563
namespace std {
template <>
struct is_constructible<std::reference_wrapper<const std::any>,
                        const std::reference_wrapper<const std::any>&>
    : public true_type {};
template <>
struct is_constructible<std::reference_wrapper<const std::any>,
                        std::reference_wrapper<const std::any>&&>
    : public true_type {};
template <>
struct is_constructible<std::reference_wrapper<const std::any>,
                        std::reference_wrapper<const std::any>&>
    : public true_type {};
template <>
struct is_copy_constructible<std::reference_wrapper<const std::any>>
    : public true_type {};
}  // namespace std

#endif

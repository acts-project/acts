// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <type_traits>

namespace Acts::detail {

namespace {
// clang-format off
    template <bool... values>
    struct all_of : std::true_type {};

    template <bool... others>
    struct all_of<false, others...> : public std::false_type {};

    template <bool... others>
    struct all_of<true, others...> : public all_of<others...> {};
// clang-format on
}  // end of anonymous namespace

template <bool... values>
constexpr bool all_of_v = all_of<values...>::value;
}  // namespace Acts::detail

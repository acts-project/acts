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
    struct any_of : std::false_type {};

    template <bool... others>
    struct any_of<true, others...> : public std::true_type {};

    template <bool... others>
    struct any_of<false, others...> : public any_of<others...> {};
// clang-format on
}  // namespace

template <bool... values>
constexpr bool any_of_v = any_of<values...>::value;
}  // namespace Acts::detail

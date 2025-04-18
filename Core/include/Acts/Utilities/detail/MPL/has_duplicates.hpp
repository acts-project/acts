// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
namespace Acts::detail {

namespace {
template <typename... Args>
struct has_duplicates;

template <>
struct has_duplicates<> {
  static constexpr bool value = false;
};

template <typename last>
struct has_duplicates<last> {
  static constexpr bool value = false;
};

template <typename first, typename second, typename... others>
struct has_duplicates<first, second, others...> {
 private:
  static constexpr bool _first = has_duplicates<first, others...>::value;
  static constexpr bool _second = has_duplicates<second, others...>::value;

 public:
  static constexpr bool value = _first || _second;
};

template <typename first, typename... others>
struct has_duplicates<first, first, others...> {
  static constexpr bool value = true;
};
}  // namespace

template <typename... Args>
constexpr bool has_duplicates_v = has_duplicates<Args...>::value;
}  // namespace Acts::detail

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/HashedString.hpp"

#include <concepts>
#include <ranges>
#include <string_view>
#include <typeindex>

namespace Acts {

/// Description of one dynamic track or track-state column.
struct DynamicColumnInfo {
  HashedString key;       ///< Key used for component access.
  std::string_view name;  ///< Original name, borrowed from the backend.
  std::type_index type;   ///< C++ value type stored in the column.
};

namespace detail {
template <typename R>
concept DynamicColumnInfoRange =
    std::ranges::input_range<R> &&
    std::same_as<std::ranges::range_value_t<R>, DynamicColumnInfo>;

inline auto dynamicColumnView(const auto& columns) {
  return columns | std::views::transform([](const auto& entry) {
           const auto& [key, column] = entry;
           return DynamicColumnInfo{key, column->name(), column->type()};
         });
}
}  // namespace detail
}  // namespace Acts

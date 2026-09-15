// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/algebra/common/ksm/concepts.hpp"

// System include(s)
#include <cstddef>
#include <tuple>

// This file describes how a single row of a substructure looks.
namespace detray::ksm {
/// @brief A single row of a substructure.
///
/// A row consists of a series of values, each of which may be an integral
/// constant, a variable, or a symbolic operation. It carries no storage and no
/// scalar type, just as a substructure does.
///
/// Note that in some cases, we might also use "row" to describe a column,
/// like when we are performing a dot product. In that case, we simply take
/// a column and transpose it into a row to avoid too much abuse of
/// terminology.
template <typename... value_ts>
struct row {
  static constexpr std::size_t num_values = sizeof...(value_ts);

  static_assert((concepts::is_value<value_ts> && ...));

  /// @brief Get the value at index I of this row.
  template <std::size_t I>
  using value_at = std::tuple_element<I, std::tuple<value_ts...>>::type;
};
}  // namespace detray::ksm

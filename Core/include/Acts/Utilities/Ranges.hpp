// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <ranges>

namespace Acts::Ranges {

/// Adaptor for converting a range to a container
/// @tparam Container the container type to convert to
template <template <typename...> class Container>
struct to_adaptor {
  /// Convert a range to a container
  /// @tparam Range the range type to convert
  /// @param range the range to convert
  /// @return the converted container
  template <typename Range>
  auto operator()(Range&& range) const {
    using ValueType = std::ranges::range_value_t<Range>;
    return Container<ValueType>(std::ranges::begin(range),
                                std::ranges::end(range));
  }
};

/// Overload operator| for piping
/// @tparam Range the range type to convert
/// @tparam Container the container type to convert to
/// @param range the range to convert
/// @param adaptor the adaptor to use
/// @return the converted container
template <typename Range, template <typename...> class Container>
auto operator|(Range&& range, to_adaptor<Container> adaptor) {
  return adaptor(std::forward<Range>(range));
}

/// Create the adaptor objects
/// @tparam Container the container type to convert to
/// @return the adaptor object
template <template <typename...> class Container>
inline constexpr to_adaptor<Container> to{};

}  // namespace Acts::Ranges

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
template <template <typename...> class Container>
struct to_adaptor {
  template <typename Range>
  auto operator()(Range&& range) const {
    using ValueType = std::ranges::range_value_t<Range>;
    return Container<ValueType>(std::ranges::begin(range),
                                std::ranges::end(range));
  }
};

// Overload operator| for piping
template <typename Range, template <typename...> class Container>
auto operator|(Range&& range, to_adaptor<Container> adaptor) {
  return adaptor(std::forward<Range>(range));
}

// Create the adaptor objects
template <template <typename...> class Container>
inline constexpr to_adaptor<Container> to{};

}  // namespace Acts::Ranges

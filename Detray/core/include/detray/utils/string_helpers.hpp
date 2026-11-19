// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/containers.hpp"

// System include(s)
#include <ranges>
#include <string>
#include <string_view>

namespace detray::utils {

/// @brief Convenience class to statically concatenate two string views.
struct string_view_concat2 {
  /// First string view in the concatenation.
  std::string_view s1;
  /// Second string view in the concatenation.
  std::string_view s2;

  /// Materialize the two views as one owning string.
  explicit operator std::string() const {
    return std::string(s1) + std::string(s2);
  }
};

/// Split string @p input at every occurrence of @p delim
/// @param input string to split
/// @param delim character that separates tokens
inline dvector<std::string> split_at_delim(const std::string &input,
                                           const char delim) {
  dvector<std::string> tokens{};

  for (const auto char_range : std::views::split(input, delim)) {
    std::string s{""};
    // TODO: Remove when range constructor becomes available in c++23
    for (const char c : char_range) {
      s.push_back(c);
    }
    tokens.push_back(std::move(s));
  }

  return tokens;
}

}  // namespace detray::utils

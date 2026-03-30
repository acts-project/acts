// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <format>
#include <ranges>

namespace Acts {

/// Utuility to join a range of strings with a delimiter.
/// Accepts any range of elements convertible to `std::string_view`.
/// @param strings Range of strings to join
/// @param delimiter Delimiter to insert between elements
/// @returns Joined string
template <std::ranges::range R>
  requires std::convertible_to<std::ranges::range_value_t<R>, std::string_view>
std::string joinStrings(R&& strings, std::string_view delimiter) {
  std::string result;
  bool first = true;

  for (auto&& item :
       strings | std::views::transform(
                     [](const auto& s) -> std::string_view { return s; })) {
    if (!first) {
      result += delimiter;
    }
    result += item;
    first = false;
  }

  return result;
}

namespace detail {
/// This mimics the signature of C++23's std::formattable concept
template <typename T, typename>
concept formattable = requires(T t) { std::format("{}", t); };
}  // namespace detail

/// Utility to join a range of formattable elements with a delimiter and custom
/// format string.
/// @param values Range of values to join
/// @param delimiter Delimiter to insert between elements
/// @param format Format string to apply to each element
/// @returns Joined string
template <std::ranges::range R>
  requires detail::formattable<std::ranges::range_value_t<R>, char>
std::string joinStrings(
    R&& values, std::string_view delimiter,
    std::format_string<const std::ranges::range_value_t<R>&> format) {
  std::string result;
  bool first = true;

  for (const auto& value : values) {
    if (!first) {
      result += delimiter;
    }
    result += std::format(format, value);
    first = false;
  }

  return result;
}

/// @cond
template <std::ranges::range R>
  requires(
      detail::formattable<std::ranges::range_value_t<R>, char> &&
      !std::convertible_to<std::ranges::range_value_t<R>, std::string_view>)
std::string joinStrings(R&& values, std::string_view delimiter) {
  std::string result;
  bool first = true;

  for (const auto& value : values) {
    if (!first) {
      result += delimiter;
    }
    result += std::format("{}", value);
    first = false;
  }

  return result;
}
/// @endcond

}  // namespace Acts

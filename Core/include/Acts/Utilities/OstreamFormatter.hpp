// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <format>
#include <sstream>
#include <string_view>

namespace Acts::detail {

template <typename Char>
struct BasicOstreamFormatter
    : std::formatter<std::basic_string_view<Char>, Char> {
  template <typename T, typename OutputIt>
  auto format(const T& value, std::basic_format_context<OutputIt, Char>& ctx)
      const -> OutputIt {
    std::basic_stringstream<Char> ss;
    ss << value;
    return std::formatter<std::basic_string_view<Char>, Char>::format(ss.view(),
                                                                      ctx);
  }
};

using OstreamFormatter = BasicOstreamFormatter<char>;

}  // namespace Acts::detail

#define ACTS_OSTREAM_FORMATTER(type) \
  template <>                        \
  struct std::formatter<type> : Acts::detail::OstreamFormatter {}

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <string>

#include <boost/algorithm/string/trim.hpp>

namespace ActsExamples::detail {

/// Return the first top-level template argument of a type name, e.g.
/// "vector<float>" -> "float" and "vector<vector<int> >" -> "vector<int>".
/// Returns an empty string if the name has no template argument.
inline std::string firstTemplateArg(const std::string& s) {
  auto lt = s.find('<');
  if (lt == std::string::npos) {
    return "";
  }
  int depth = 0;
  std::size_t start = lt + 1;
  for (std::size_t i = lt; i < s.size(); ++i) {
    char c = s[i];
    if (c == '<') {
      ++depth;
      if (depth == 1) {
        start = i + 1;
      }
    } else if (c == '>') {
      --depth;
      if (depth == 0) {
        return boost::algorithm::trim_copy(s.substr(start, i - start));
      }
    } else if (c == ',' && depth == 1) {
      return boost::algorithm::trim_copy(s.substr(start, i - start));
    }
  }
  return "";
}

}  // namespace ActsExamples::detail

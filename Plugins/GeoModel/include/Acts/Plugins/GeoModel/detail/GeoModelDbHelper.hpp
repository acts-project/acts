// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <string>
#include <vector>

namespace Acts::detail::GeoModelDbHelper {

/// @brief Split a string into a vector of strings
///
/// @param entry the string to be split
/// @param deliminater split indicator
///
/// @return a vector of strings
inline std::vector<std::string> tokenize(const std::string& entry,
                                         const std::string& deliminater) {
  std::vector<std::string> result;
  std::string currentEntry = entry;
  std::size_t pos = 0;
  while ((pos = currentEntry.find(deliminater)) != std::string::npos) {
    auto found = currentEntry.substr(0, pos);
    if (!found.empty()) {
      result.push_back(currentEntry.substr(0, pos));
    }
    currentEntry.erase(0, pos + deliminater.length());
  }
  if (!currentEntry.empty()) {
    result.push_back(currentEntry);
  }
  return result;
}

}  // namespace Acts::detail::GeoModelDbHelper

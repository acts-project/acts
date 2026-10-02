// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <optional>
#include <set>
#include <string>
#include <vector>

#include <traccc/definitions/primitives.hpp>
#include <traccc/edm/particle.hpp>
#include <traccc/efficiency/track_matcher.hpp>

namespace traccc {
stepped_percentage::stepped_percentage(scalar ratio) : m_min_ratio(ratio) {}

std::string stepped_percentage::get_name() const {
  char buffer[512];
  snprintf(buffer, 512, "Stepped with ≥ %.1f%% similarity",
           100.f * m_min_ratio);
  return std::string(buffer);
}

std::optional<std::uint64_t> stepped_percentage::operator()(
    const std::vector<std::vector<std::uint64_t>>& p) const {
  std::multiset<std::uint64_t> cnt;

  /*
   * Record all the particle identifiers in a multiset in order to count
   * them.
   */
  for (const std::vector<std::uint64_t>& i : p) {
    for (std::uint64_t j : i) {
      cnt.insert(j);
    }
  }

  /*
   * From the maximum size, decrease the matching count required until we
   * find a match or we dip below the threshold.
   */
  for (std::size_t n = p.size();
       n <= p.size() &&
       (static_cast<float>(n) / static_cast<float>(p.size())) > m_min_ratio;
       --n) {
    for (std::uint64_t i : cnt) {
      if (cnt.count(i) == n) {
        return {i};
      }
    }
  }

  return {};
}

exact::exact() {}

std::string exact::get_name() const {
  char buffer[512];
  snprintf(buffer, 512, "Exact");
  return std::string(buffer);
}

std::optional<std::uint64_t> exact::operator()(
    const std::vector<std::vector<std::uint64_t>>& p) const {
  std::multiset<std::uint64_t> cnt;

  /*
   * Record all the particle identifiers in a multiset in order to count
   * them.
   */
  for (const std::vector<std::uint64_t>& i : p) {
    for (std::uint64_t j : i) {
      cnt.insert(j);
    }
  }

  /*
   * Find a particle which matches every single spacepoint.
   */
  for (std::uint64_t i : cnt) {
    if (cnt.count(i) == p.size()) {
      return {i};
    }
  }

  return {};
}
}  // namespace traccc

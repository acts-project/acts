// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"

// System include(s)
#include <string>
#include <string_view>
#include <unordered_map>

namespace detray {

/// Map volume indices to volume names and vice versa
struct name_map {
  /// @brief set the name of the detector
  DETRAY_HOST void set_detector_name(const std::string_view name) {
    detector_name = name;
  }

  /// @returns the name of the detector
  DETRAY_HOST const std::string& get_detector_name() const {
    return detector_name;
  }

  /// @returns whether the name map is empty
  DETRAY_HOST bool empty() const {
    assert(index_to_name.size() == name_to_index.size());
    return index_to_name.empty();
  }

  /// @returns whether the name map contains a given volume index
  DETRAY_HOST bool contains(dindex idx) const {
    return index_to_name.contains(idx);
  }

  /// @returns whether the name map contains a given volume index
  DETRAY_HOST bool contains(std::string_view name) const {
    return name_to_index.contains(name);
  }

  /// @brief emplace a new index <-> name mapping
  DETRAY_HOST void emplace(dindex idx, const std::string& name) {
    index_to_name.try_emplace(idx, name);
    name_to_index.try_emplace(index_to_name.at(idx), idx);
  }

  /// @returns the name of a volume given the volume index @param idx - const
  DETRAY_HOST decltype(auto) at(dindex idx) const {
    return index_to_name.at(idx);
  }

  /// @returns the name of a volume given the volume index @param idx
  DETRAY_HOST decltype(auto) at(dindex idx) { return index_to_name.at(idx); }

  /// @returns the index of a volume given the volume name @param name
  DETRAY_HOST decltype(auto) at(const std::string& name) const {
    return name_to_index.at(name);
  }

  /// @returns the index of a volume given the volume name @param name
  DETRAY_HOST decltype(auto) at(const std::string_view name) const {
    return name_to_index.at(name);
  }

  /// Clears the name map
  DETRAY_HOST void clear() {
    detector_name = "";
    index_to_name.clear();
    name_to_index.clear();
  }

  /// Clears the volume name mapping, but keeps detector name
  DETRAY_HOST void clear_names() {
    index_to_name.clear();
    name_to_index.clear();
  }

 private:
  /// Name of the detector
  std::string detector_name{};
  /// Map the volume names to volume indices (position of the volume
  /// descriptors in the detector volume container)
  /// @ {
  std::unordered_map<dindex, std::string> index_to_name{};
  std::unordered_map<std::string_view, dindex> name_to_index{};
  /// @}
};

}  // namespace detray

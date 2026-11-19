// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include "detray/detectors/metadata_registry.hpp"

namespace detray::registry {

/// Stable family ID for a registered metadata type; empty for unregistered
/// types.
template <concepts::metadata metadata_t>
constexpr std::string_view metadata_id() {
  using algebra_t = typename metadata_t::algebra_type;
  // clang-format off
#define DETRAY_MATCH_METADATA(ID) \
  if constexpr (std::is_same_v<metadata_t, ID##_metadata<algebra_t>>) { return #ID; }
  DETRAY_FOR_EACH_METADATA(DETRAY_MATCH_METADATA)
#undef DETRAY_MATCH_METADATA
  // clang-format on
  return {};
}

/// Accept stable IDs and the previous '<id>_metadata' spelling from JSON files.
inline std::string_view resolve_metadata_id(std::string_view id) {
  for (auto known : metadata_ids) {
    if (id == known) {
      return known;
    }
  }
  constexpr std::string_view suffix = "_metadata";
  if (id.ends_with(suffix)) {
    id.remove_suffix(suffix.size());
  }
  for (auto known : metadata_ids) {
    if (id == known) {
      return known;
    }
  }
  throw std::invalid_argument("Unknown or unavailable detector metadata ID '" +
                              std::string(id) +
                              "'; supply an available metadata ID explicitly");
}

namespace detail {
template <typename algebra_t, std::size_t I = 0, typename callable_t>
decltype(auto) visit_metadata_impl(std::string_view id, callable_t&& callable) {
  using metadata_t = std::tuple_element_t<I, metadata_types<algebra_t>>;
  if (id == metadata_ids[I]) {
    return std::forward<callable_t>(callable).template operator()<metadata_t>();
  }
  if constexpr (I + 1 < metadata_ids.size()) {
    return visit_metadata_impl<algebra_t, I + 1>(
        id, std::forward<callable_t>(callable));
  } else {
    throw std::invalid_argument("Unavailable detector metadata ID: " +
                                std::string(id));
  }
}
}  // namespace detail

/// Dispatch a metadata ID to a concrete algebra-parametric family.
template <typename algebra_t, typename callable_t>
decltype(auto) visit_metadata(std::string_view id, callable_t&& callable) {
  return detail::visit_metadata_impl<algebra_t>(
      resolve_metadata_id(id), std::forward<callable_t>(callable));
}

}  // namespace detray::registry

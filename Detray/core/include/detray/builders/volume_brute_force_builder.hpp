// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/builders/volume_accelerator_builder_interface.hpp"
#include "detray/core/concepts.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/utils/logging.hpp"

// System include(s)
#include <algorithm>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace detray {

/// @brief Builds volume brute force accelerators for the daughter volumes of
/// root volumes.
///
/// Every registered root gets one entry in the volume brute force collection
/// (@c accel_id::e_volume_brute_force ), which contains the indices of its
/// daughter volumes. The root links to the entry through its @c e_volume
/// accelerator link.
///
/// If the volume brute force accelerator is also the default volume
/// accelerator, the builder reserves index 0 for the world volume lookup: It
/// contains the daughters of the root @c dindex_invalid , or all top-level
/// volumes, if this root was not registered.
///
/// The builder only checks properties that are specific to the brute force
/// accelerator. The volume hierarchy is checked by the consistency checker.
///
/// @tparam detector_t the type of the detector
/// @note The template parameter is not constrained, because the detector
/// declares the builder as friend while the detector type is incomplete.
template <typename detector_t>
class volume_brute_force_builder final
    : public volume_accelerator_builder_interface<detector_t> {
  static_assert(concepts::has_volume_brute_force<detector_t>,
                "Detector type does not contain a volume brute force "
                "accelerator");

  /// Accelerator ids of the detector
  using accel_id = typename detector_t::accel::id;
  /// The id of the volume brute force accelerator in the accelerator store
  static constexpr auto bf_id{accel_id::e_volume_brute_force};
  /// Whether the brute force accelerator holds the world volume lookup
  static constexpr bool has_world_lookup{bf_id == accel_id::e_volume_default};

 public:
  /// @returns the name of the builder
  DETRAY_HOST auto name() const -> std::string_view override {
    return "volume_brute_force_builder";
  }

  /// @copydoc volume_accelerator_builder_interface::add_root
  ///
  /// @throws std::invalid_argument if the root is already registered, if the
  /// daughter list is empty or contains duplicates, or if the root is the
  /// world, but the brute force accelerator is not the default volume
  /// accelerator
  DETRAY_HOST void add_root(const dindex root,
                            const std::vector<dindex>& daughters) override {
    std::stringstream err_stream{};
    err_stream << "Volume brute force builder: Root " << root << ": ";

    if (root == dindex_invalid && !has_world_lookup) {
      err_stream << "World volume lookup needs the default volume accelerator";
      throw_error(err_stream);
    }
    if (m_roots.contains(root)) {
      err_stream << "Root registered twice";
      throw_error(err_stream);
    }
    if (daughters.empty()) {
      err_stream << "No daughter volumes";
      throw_error(err_stream);
    }
    std::vector<dindex> sorted{daughters};
    std::ranges::sort(sorted);
    if (std::ranges::adjacent_find(sorted) != sorted.end()) {
      err_stream << "Duplicate daughter volumes";
      throw_error(err_stream);
    }

    m_roots.emplace(root, daughters);
  }

  /// @copydoc volume_accelerator_builder_interface::roots
  DETRAY_HOST auto roots() const -> std::vector<dindex> override {
    std::vector<dindex> roots{};
    roots.reserve(m_roots.size());
    for (const auto& [root, _] : m_roots) {
      roots.push_back(root);
    }
    return roots;
  }

  /// @copydoc volume_accelerator_builder_interface::build
  ///
  /// @throws std::invalid_argument if a root does not exist in the detector,
  /// if its daughter volumes are already registered, or if the world root is
  /// registered, but another builder already built the world volume lookup
  DETRAY_HOST void build(detector_t& det) override {
    constexpr auto vol_obj_id{detector_t::geo_obj_ids::e_volume};
    auto& bf_coll = det._accelerators.template get<bf_id>();

    // Index 0 of the default volume accelerator: World volume lookup
    if constexpr (has_world_lookup) {
      const auto search = m_roots.find(dindex_invalid);
      if (search != m_roots.end() && bf_coll.size() != 0u) {
        std::stringstream err_stream{};
        err_stream << "Volume brute force builder: World volume lookup was "
                      "already built by another builder: Register the world "
                      "root with the first volume brute force builder";
        throw_error(err_stream);
      }
      if (bf_coll.size() == 0u) {
        if (search != m_roots.end()) {
          bf_coll.push_back(search->second);
        } else {
          std::vector<dindex> top_level{};
          for (const auto& vol_desc : det.volumes()) {
            if (!vol_desc.has_mother()) {
              top_level.push_back(vol_desc.index());
            }
          }
          bf_coll.push_back(top_level);
        }
      }
    }

    for (const auto& [root, daughters] : m_roots) {
      if (root == dindex_invalid) {
        continue;
      }

      std::stringstream err_stream{};
      err_stream << "Volume brute force builder: Root " << root << ": ";
      if (root >= det.volumes().size()) {
        err_stream << "Volume does not exist in detector";
        throw_error(err_stream);
      }
      auto& root_desc = det._volumes.at(root);
      if (!root_desc.template accel_link<vol_obj_id>().is_invalid()) {
        err_stream << "Daughter volumes are already registered";
        throw_error(err_stream);
      }

      bf_coll.push_back(daughters);
      root_desc.template set_accel_link<vol_obj_id>(
          bf_id, static_cast<dindex>(bf_coll.size() - 1u));

      DETRAY_VERBOSE_HOST("-> Volume " << root << ": " << daughters.size()
                                       << " daughter volumes");
    }
  }

 private:
  /// Log and throw an error
  /// @param err_stream the error message
  DETRAY_HOST static void throw_error(const std::stringstream& err_stream) {
    DETRAY_FATAL_HOST(err_stream.str());
    throw std::invalid_argument(err_stream.str());
  }

  /// Daughter volumes per root volume
  std::map<dindex, std::vector<dindex>> m_roots{};
};

}  // namespace detray

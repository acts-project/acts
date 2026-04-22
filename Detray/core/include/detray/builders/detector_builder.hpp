// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/builders/volume_builder.hpp"
#include "detray/builders/volume_builder_interface.hpp"
#include "detray/core/concepts.hpp"
#include "detray/core/detector.hpp"
#include "detray/utils/detector_statistics.hpp"
#include "detray/utils/logging.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace detray {

/// @brief Provides functionality to build a detray detector volume by volume
///
/// @tparam metadata the type definitions for the detector
/// @tparam volume_builder_t the basic volume builder to be used for the
///                          geometry data
/// @tparam volume_data_t the data structure that holds the volume builders
template <typename metadata,
          template <typename> class volume_builder_t = volume_builder,
          template <typename...> class volume_data_t = std::vector>
class detector_builder {
 public:
  using detector_type = detector<metadata, host_container_types>;
  using algebra_type = typename detector_type::algebra_type;
  using scalar_type = dscalar<algebra_type>;

  DETRAY_HOST detector_builder() { DETRAY_VERBOSE_HOST("New builder created"); }

  /// Set the name of the detector under construction to @param det_name
  DETRAY_HOST void set_name(std::string det_name) {
    m_detector_name = std::move(det_name);
    DETRAY_VERBOSE_HOST("Set detector name: " << m_detector_name);
  }

  /// @returns the name of the detector under construction
  DETRAY_HOST std::string_view name() const { return m_detector_name; }

  /// Add a new volume builder that will build a volume of the shape given by
  /// @param id
  template <typename... Args>
  DETRAY_HOST auto new_volume(const volume_id id, Args&&... args)
      -> volume_builder_interface<detector_type>* {
    DETRAY_VERBOSE_HOST("Adding new volume: " << m_volumes.size());

    m_volumes.push_back(std::make_unique<volume_builder_t<detector_type>>(
        id, static_cast<dindex>(m_volumes.size()),
        std::forward<Args>(args)...));

    return m_volumes.back().get();
  }

  /// @returns the number of volumes currently registered in the builder
  DETRAY_HOST auto n_volumes() const -> dindex {
    return static_cast<dindex>(m_volumes.size());
  }

  /// @returns 'true' if there is a volume builder registered for
  /// the volume with index @param volume_idx
  DETRAY_HOST bool has_volume(const std::size_t volume_idx) const {
    return volume_idx < m_volumes.size();
  }

  /// Decorate a volume builder at position @param volume_idx with more
  /// functionality
  template <class builder_t>
  DETRAY_HOST auto decorate(dindex volume_idx) -> builder_t* {
    assert(has_volume(volume_idx));

    m_volumes[volume_idx] =
        std::make_unique<builder_t>(std::move(m_volumes[volume_idx]));

    // Always works, we set it as this type in the line above
    return dynamic_cast<builder_t*>(m_volumes[volume_idx].get());
  }

  /// Decorate a volume builder @param v_builder with more functionality
  template <class builder_t>
  DETRAY_HOST auto decorate(
      const volume_builder_interface<detector_type>* v_builder) -> builder_t* {
    assert(v_builder != nullptr);

    return decorate<builder_t>(v_builder->vol_index());
  }

  /// Access a particular volume builder by volume index @param volume_idx
  DETRAY_HOST
  auto operator[](dindex volume_idx)
      -> volume_builder_interface<detector_type>* {
    return m_volumes[volume_idx].get();
  }

  /// Assembles the final detector from the volumes builders and allocates
  /// the detector containers with the memory resource @param resource
  DETRAY_HOST
  auto build(vecmem::memory_resource& resource) -> detector_type {
    DETRAY_INFO_HOST("Building detector: \"" << name() << "\"... ");
    DETRAY_INFO_HOST("-> type: " << DETRAY_TYPENAME(metadata));

    detector_type det{resource};

    DETRAY_VERBOSE_HOST("Have " << m_volumes.size()
                                << " configured volume builders");
    DETRAY_VERBOSE_HOST("Start building the volumes...");
    for (auto& vol_builder : m_volumes) {
      DETRAY_VERBOSE_HOST("-> Build: " << vol_builder->name());
      vol_builder->build(det);
    }

    // TODO: Make fully generic for more volume accelerator types

    // TODO: Add sorting, data deduplication etc. here later...

    DETRAY_INFO_HOST("-> Built " << det.volumes().size() << " volumes");
    DETRAY_INFO_HOST("-> Built " << det.surfaces().size() << " surfaces:");
    DETRAY_INFO_HOST("--> portals:    " << detray::n_portals(det));
    DETRAY_INFO_HOST("--> sensitives: " << detray::n_sensitives(det));
    DETRAY_INFO_HOST("--> passives:   " << detray::n_passives(det));

    if constexpr (detray::concepts::has_surface_grids<detector_type>) {
      DETRAY_INFO_HOST("-> Built " << detray::n_surface_grids(det)
                                   << " surface grids");
    }

    if constexpr (detray::concepts::has_material_maps<detector_type>) {
      DETRAY_INFO_HOST("-> Built " << detray::n_material_maps(det)
                                   << " material maps");
    }

    if constexpr (detray::concepts::has_homogeneous_material<detector_type>) {
      DETRAY_INFO_HOST("-> Built homogeneous material:");
      DETRAY_INFO_HOST("--> slabs: " << detray::n_material_slabs(det));
      DETRAY_INFO_HOST("--> rods:  " << detray::n_material_rods(det));
    }
    DETRAY_INFO_HOST("Detector building complete: " << name());

    return det;
  }

  /// Assembles the final detector and fill an externally provided name map
  /// @param name_map
  DETRAY_HOST
  auto build(vecmem::memory_resource& resource,
             typename detector_type::name_map& name_map) -> detector_type {
    DETRAY_VERBOSE_HOST("detray: filling names for detector " << name());

    assert(name_map.empty());

    // By convention the name of the detector is at position 0
    name_map.set_detector_name(m_detector_name);

    for (auto& vol_builder : m_volumes) {
      name_map.emplace(vol_builder->vol_index(),
                       std::string{vol_builder->name()});
    }

    return build(resource);
  }

 private:
  /// Name of the new detector
  std::string m_detector_name{"detray_detector"};
  /// Data structure that holds a volume builder for every detector volume
  volume_data_t<std::unique_ptr<volume_builder_interface<detector_type>>>
      m_volumes{};
};

}  // namespace detray

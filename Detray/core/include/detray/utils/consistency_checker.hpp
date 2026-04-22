// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/geometry/tracking_volume.hpp"
#include "detray/material/concepts.hpp"
#include "detray/material/predefined_materials.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/ranges.hpp"
#include "detray/utils/type_registry.hpp"

// System include(s)
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace detray::detail {

/// Checks every collection in a multistore to be empty and prints a warning
template <typename store_t, std::size_t... I>
void report_empty(const store_t &store,
                  [[maybe_unused]] const std::string &store_name,
                  std::index_sequence<I...> /*seq*/) {
  ((store.template empty<types::id_cast<typename store_t::value_types, I>>()
        ?
#if DETRAY_LOG_LVL < 0
        std::clog << ""
#else
        // The log macro does not compile here...
        std::clog << "DETRAY WARNING (HOST): " << __FILENAME__ << ":"
                  << __LINE__ << " " << store_name
                  << " has empty collection no. " << I << std::endl
#endif
        : std::clog << ""),
   ...);
}

/// A functor that checks the surface descriptor and correct volume index in
/// every acceleration data structure for a given volume
struct surface_checker {
  /// Test the contained surfaces for consistency
  template <typename detector_t>
  DETRAY_HOST_DEVICE void operator()(
      const typename detector_t::surface_type &sf_descr, const detector_t &det,
      const dindex vol_idx, const typename detector_t::name_map &names) const {
    const auto sf = geometry::surface{det, sf_descr};
    const auto vol = tracking_volume{det, vol_idx};
    std::stringstream err_stream{};
    err_stream << "Volume \"" << vol.name(names) << "\":\n";

    if (!sf.self_check(err_stream)) {
      throw std::invalid_argument(err_stream.str());
    }

    if (sf.volume() != vol_idx) {
      err_stream << "Incorrect volume index on surface: vol. index " << vol_idx
                 << ", sf: " << sf;

      DETRAY_FATAL_HOST(err_stream.str());
      throw std::invalid_argument(err_stream.str());
    }

    // Does the mask link to an existing volume?
    const auto vol_links = sf.volume_links();
    for (const auto vol_link : vol_links) {
      if (!detail::is_invalid_value(vol_link) &&
          (vol_link >= det.volumes().size())) {
        err_stream << "Incorrect volume link to non-existent volume "
                   << vol_link << " on surface " << sf;

        DETRAY_FATAL_HOST(err_stream.str());
        throw std::invalid_argument(err_stream.str());
      }
    }

    // A passive surface should have material, if the detector was
    // constructed with material
    if (!det.material_store().all_empty() &&
        (sf.is_passive() && !sf.has_material())) {
      DETRAY_WARN_HOST(err_stream.str());
      DETRAY_WARN_HOST("Passive surface without material: " << sf_descr
                                                            << "\n");
    }

    // Check that the same surface is registered in the detector surface
    // lookup
    const auto sf_from_lkp =
        geometry::surface{det, det.surface(sf.identifier())};
    if (sf_from_lkp != sf) {
      err_stream << "Surfaces in volume and detector lookups "
                 << "differ:\n In volume: " << vol
                 << "\nFound in detector surface lookup: " << sf_from_lkp;

      DETRAY_FATAL_HOST(err_stream.str());
      throw std::runtime_error(err_stream.str());
    }
  }

  /// Test whether a given surface @param check_desc is properly registered at
  /// least once in one of the volume acceleration data structures
  ///
  /// @param ref_descr one of the surfaces in the volumes acceleration data
  /// @param check_descr surface that we are searching for
  /// @param success communicate success to the outside
  template <typename sf_descr_t, typename sf_source_descr_t>
  DETRAY_HOST_DEVICE void operator()(const sf_descr_t &ref_descr,
                                     const sf_source_descr_t &check_descr,
                                     bool &success) const {
    // Check that the other surfaces in the acceleration structure belong
    // there The volume index of the check_descr must be checked to be
    // correct beforehand, e.g. by the call operator above
    if (ref_descr.volume() != check_descr.volume()) {
      std::stringstream err_stream{};
      err_stream << "Inconsistent volume index in accel. structure: "
                 << "Expected: " << check_descr.volume()
                 << ". Got surface in accel. structure: " << ref_descr;

      DETRAY_FATAL_HOST(err_stream.str());
      throw std::invalid_argument(err_stream.str());
    }

    // Check if it is the surface we are looking for
    if (ref_descr == check_descr) {
      success = true;
    }
  }
};

/// A functor that checks the material parametrization for a surface/volume
struct material_checker {
  /// Error message for material consistency check
  template <typename material_t>
  void throw_material_error(const std::string &type, const dindex idx,
                            const material_t &mat) const {
    std::stringstream err_stream{};
    err_stream << "Invalid material found in: " << type << " at index " << idx
               << ": " << mat;
    DETRAY_FATAL_HOST(err_stream.str());
    throw std::invalid_argument(err_stream.str());
  }

  /// Test whether a given material map contains invalid material
  ///
  /// @param material_coll collection of material grids
  /// @param idx the specific grid to be checked
  /// @param id type id of the material grid collection
  template <typename material_coll_t, typename index_t, typename id_t>
    requires concepts::material_map<typename material_coll_t::value_type>
  DETRAY_HOST_DEVICE void operator()(const material_coll_t &material_coll,
                                     const index_t idx, const id_t id) const {
    try {
      const auto mat_map = material_coll.at(idx);

      // Check whether there are any entries in the bins
      if (mat_map.size() == 0u) {
        std::stringstream err_stream{};
        err_stream << "Empty material grid: " << static_cast<int>(id)
                   << " at index " << idx;
        DETRAY_FATAL_HOST(err_stream.str());
        throw std::invalid_argument(err_stream.str());
      } else {
        for (const auto &bin : mat_map.bins()) {
          if (bin.size() == 0u) {
            std::stringstream err_stream{};
            err_stream << "Empty material bin: " << static_cast<int>(id)
                       << " at index " << idx;
            DETRAY_FATAL_HOST(err_stream.str());
            throw std::invalid_argument(err_stream.str());
          }
        }
      }
    } catch (std::out_of_range &) {
      std::stringstream err_stream{};
      err_stream << "Out of range material access in: "
                 << "binned material collection at index " << idx;
      DETRAY_FATAL_HOST(err_stream.str());
      throw std::invalid_argument(err_stream.str());
    }
  }

  /// Test whether a given collection of material contains invalid material
  ///
  /// @param material_coll collection of material slabs/rods/raw mat
  /// @param idx the specific instance to be checked
  template <typename material_coll_t, typename index_t, typename id_t>
    requires concepts::homogeneous_material<
        typename material_coll_t::value_type>
  DETRAY_HOST_DEVICE void operator()(const material_coll_t &material_coll,
                                     const index_t idx, const id_t) const {
    using material_t = typename material_coll_t::value_type;
    using scalar_t = typename material_t::scalar_type;

    try {
      const material_t &mat = material_coll.at(idx);

      // Homogeneous volume material
      if constexpr (concepts::material_params<material_t>) {
        if (mat == detray::vacuum<scalar_t>{}) {
          throw_material_error("homogeneous volume material", idx, mat);
        }

      } else {
        // Material slabs and rods
        if (!mat) {
          throw_material_error("homogeneous surface material", idx, mat);
        }
      }
    } catch (std::out_of_range &) {
      std::stringstream err_stream{};
      err_stream << "Out of range material access in: "
                 << "homogeneous material collection at index " << idx;
      DETRAY_FATAL_HOST(err_stream.str());
      throw std::invalid_argument(err_stream.str());
    }
  }
};

/// @brief Checks whether the data containers of a detector are empty
///
/// In case the default metadata is used, the unused containers are allowed to
/// be empty.
template <typename detector_t>
inline void check_empty(const detector_t &det, const bool verbose) {
  // Check if there is at least one portal in the detector
  auto find_portals = [&det]() {
    if (det.portals().empty()) {
      return false;
    }
    // In the brute force finder, also other surfaces can be contained, e.g.
    // passive surfaces (depends on the detector)
    return std::ranges::any_of(
        det.portals(), [](auto pt_desc) { return pt_desc.is_portal(); });
  };

  // TODO: Check for empty volume acceleration structures

  // Fatal errors
  if (det.volumes().empty()) {
    std::string err_str{"No volumes in detector"};
    DETRAY_FATAL_HOST(err_str);
    throw std::runtime_error(err_str);
  }
  if (det.surfaces().empty()) {
    std::string err_str{"No surfaces found"};
    DETRAY_FATAL_HOST(err_str);
    throw std::runtime_error(err_str);
  }
  if (det.transform_store().empty()) {
    std::string err_str{"No transforms in detector"};
    DETRAY_FATAL_HOST(err_str);
    throw std::runtime_error(err_str);
  }
  if (det.mask_store().all_empty()) {
    std::string err_str{"No masks in detector"};
    DETRAY_FATAL_HOST(err_str);
    throw std::runtime_error(err_str);
  }
  if (!find_portals()) {
    std::string err_str{"No portals in detector"};
    DETRAY_FATAL_HOST(err_str);
    throw std::runtime_error(err_str);
  }

  // Warnings

  // Check the material description
  if (det.material_store().all_empty()) {
    DETRAY_WARN_HOST("No material in detector");
  } else if (verbose) {
    // Check for empty material collections
    detail::report_empty(
        det.material_store(), "material store",
        std::make_index_sequence<detector_t::material::n_types>{});
  }

  // Check for empty acceleration data structure collections (e.g. grids)
  if (verbose) {
    // Check for empty mask collections
    detail::report_empty(
        det.mask_store(), "mask store",
        std::make_index_sequence<detector_t::masks::n_types>{});

    detail::report_empty(
        det.accelerator_store(), "acceleration data structures store",
        std::make_index_sequence<detector_t::accel::n_types>{});
  }

  // TODO: Implement volume acceleration structure check
}

/// @brief Checks the internal consistency of a detector
template <typename detector_t>
inline bool check_consistency(const detector_t &det, const bool verbose = false,
                              const typename detector_t::name_map &names = {}) {
  DETRAY_INFO_HOST("Checking detector consistency...");

  check_empty(det, verbose);

  // Check the volumes
  constexpr bool portal_eq_passive{
      detector_t::volume_type::object_id::e_portal ==
      detector_t::volume_type::object_id::e_passive};

  for (const auto &[idx, vol_desc] : detray::views::enumerate(det.volumes())) {
    const auto vol = tracking_volume{det, vol_desc};

    std::stringstream err_stream{};
    err_stream << "Volume \"" << vol.name(names) << "\":\n";

    // Check that nothing is obviously broken
    if (!vol.self_check(err_stream)) {
      DETRAY_FATAL_HOST(err_stream.str());
      throw std::invalid_argument(err_stream.str());
    }

    // Check consistency in the context of the owning detector
    if (vol.index() != idx) {
      err_stream << "Incorrect volume index! Found volume:\n"
                 << vol << "\nat index " << idx;

      DETRAY_FATAL_HOST(err_stream.str());
      throw std::invalid_argument(err_stream.str());
    }

    // Check that only surfaces belonging to this volume are in the surface
    // range of the detector

    // Check volume indices
    for (const auto sf_desc : vol.template surfaces<surface_id::e_all>()) {
      if (sf_desc.volume() != vol.index()) {
        err_stream << "Surface of other volume detected:\n\nvolume " << vol
                   << "\n\nsurface " << sf_desc;

        DETRAY_FATAL_HOST(err_stream.str());
        throw std::invalid_argument(err_stream.str());
      }
    }
    // Check surface type: portal or passive
    for (const auto sf_desc : vol.portals()) {
      if (!portal_eq_passive && !sf_desc.is_portal()) {
        err_stream << "Non-portal surface discovered in portal "
                      "range: "
                   << vol_desc.template sf_link<surface_id::e_portal>()
                   << "\n\nvolume " << vol << "\n\nsurface " << sf_desc;

        DETRAY_FATAL_HOST(err_stream.str());
        throw std::invalid_argument(err_stream.str());
      }

      if (portal_eq_passive && !(sf_desc.is_portal() || sf_desc.is_passive())) {
        err_stream << "Sensitive surface discovered in portal/passive "
                      "range: "
                   << vol_desc.template sf_link<surface_id::e_portal>()
                   << "\n\nvolume " << vol << "\n\nsurface " << sf_desc;

        DETRAY_FATAL_HOST(err_stream.str());
        throw std::invalid_argument(err_stream.str());
      }
    }
    // Check surface type: passive or portal
    for (const auto sf_desc : vol.template surfaces<surface_id::e_passive>()) {
      if (!portal_eq_passive && !sf_desc.is_passive()) {
        err_stream << "Non-passive surface discovered in passive "
                      "range: "
                   << vol_desc.template sf_link<surface_id::e_passive>()
                   << "\n\nvolume " << vol << "\n\nsurface " << sf_desc;

        DETRAY_FATAL_HOST(err_stream.str());
        throw std::invalid_argument(err_stream.str());
      }

      if (portal_eq_passive && !(sf_desc.is_portal() || sf_desc.is_passive())) {
        err_stream << "Sensitive surface discovered in portal/passive "
                      "range: "
                   << vol_desc.template sf_link<surface_id::e_passive>()
                   << "\n\nvolume " << vol << "\n\nsurface " << sf_desc;

        DETRAY_FATAL_HOST(err_stream.str());
        throw std::invalid_argument(err_stream.str());
      }
    }
    // Check surface type: sensitive
    for (const auto sf_desc :
         vol.template surfaces<surface_id::e_sensitive>()) {
      if (!sf_desc.is_sensitive()) {
        err_stream << "Non-sensitive surface discovered in sensitive "
                      "range: "
                   << vol_desc.template sf_link<surface_id::e_sensitive>()
                   << "\n\nvolume " << vol << "\n\nsurface " << sf_desc;

        DETRAY_FATAL_HOST(err_stream.str());
        throw std::invalid_argument(err_stream.str());
      }
    }

    // Go through the acceleration data structures and check the surfaces
    vol.template visit_surfaces<surface_id::e_all, detail::surface_checker>(
        det, vol.index(), names);

    // Check the volume material, if present
    if (vol.has_material()) {
      vol.template visit_material<detail::material_checker>(
          vol_desc.material().id());
    }
  }

  // Check the surfaces in the detector's surface lookup
  for (const auto &[idx, sf_desc] : detray::views::enumerate(det.surfaces())) {
    const auto sf = geometry::surface{det, sf_desc};
    const auto vol = tracking_volume{det, sf.volume()};

    std::stringstream err_stream{};
    err_stream << "Volume \"" << vol.name(names) << "\":\n";

    // Check that nothing is obviously broken
    if (!sf.self_check(err_stream)) {
      err_stream << "\nat surface no. " << std::to_string(idx);
      throw std::invalid_argument(err_stream.str());
    }

    // Check consistency in the context of the owning detector
    if (sf.index() != idx) {
      err_stream << "Incorrect surface index! Found surface:\n"
                 << sf << "\nat index " << idx;

      DETRAY_FATAL_HOST(err_stream.str());
      throw std::invalid_argument(err_stream.str());
    }

    // Check that the surface can be found in its volume's acceleration
    // data structures (if there are no grids, must at least be in the
    // brute force method)
    bool is_registered = false;
    vol.template visit_surfaces<surface_id::e_all, detail::surface_checker>(
        sf_desc, is_registered);

    if (!is_registered) {
      err_stream << "Found surface that is not part of its "
                 << "volume's navigation acceleration data structures:\n"
                 << "Surface: " << sf;

      DETRAY_FATAL_HOST(err_stream.str());
      throw std::invalid_argument(err_stream.str());
    }

    // Check the surface material, if present
    if (sf.has_material()) {
      DETRAY_DEBUG_HOST("Checking surface sf=" << sf);
      sf.template visit_material<detail::material_checker>(
          sf_desc.material().id());
    }
  }

  DETRAY_INFO_HOST("Detector consistency check: OK");

  return true;
}

}  // namespace detray::detail

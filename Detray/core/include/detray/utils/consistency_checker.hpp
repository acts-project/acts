// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/core/concepts.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/detail/vertexer.hpp"
#include "detray/geometry/shapes.hpp"
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
  /// @param sf_descr Surface descriptor being checked.
  /// @param det Detector that owns the surface.
  /// @param vol_idx Index of the volume expected to contain the surface.
  /// @param names Map used to resolve volume names in diagnostics.
  template <concepts::detector detector_t>
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

  /// Test whether a given surface @p check_desc is properly registered at
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

/// A functor that checks the daughter volumes in the volume accelerator of a
/// mother volume
struct daughter_checker {
  /// Check a daughter volume of a mother volume
  ///
  /// @tparam detector_t the type of the detector
  ///
  /// @param daughter_idx index of the daughter volume in the accelerator
  /// @param det the detector
  /// @param mother_idx index of the mother volume
  /// @param query_idx index of a volume that is searched for
  /// @param n_found number of times the queried volume was found (output)
  template <concepts::detector detector_t>
  DETRAY_HOST void operator()(const dindex daughter_idx, const detector_t &det,
                              const dindex mother_idx, const dindex query_idx,
                              dindex &n_found) const {
    std::stringstream err_stream{};
    if (daughter_idx >= det.volumes().size()) {
      err_stream << "Daughter volume index " << daughter_idx
                 << " out of bounds in mother volume " << mother_idx;
      DETRAY_FATAL_HOST(err_stream.str());
      throw std::invalid_argument(err_stream.str());
    }
    if (const dindex mother{det.volume(daughter_idx).mother()};
        mother != mother_idx) {
      err_stream << "Daughter volume " << daughter_idx
                 << " in volume accelerator of volume " << mother_idx
                 << " has a different mother volume: " << mother;
      DETRAY_FATAL_HOST(err_stream.str());
      throw std::invalid_argument(err_stream.str());
    }
    if (daughter_idx == query_idx) {
      ++n_found;
    }
  }
};

/// A functor that checks that the accelerator link of a root volume points to
/// an existing volume acceleration structure
struct volume_accel_link_checker {
  /// @tparam accel_coll_t the type of the accelerator collection
  /// @tparam accel_index_t the type of the index into the collection
  ///
  /// @param coll the collection of acceleration structures
  /// @param index the index of the acceleration structure
  ///
  /// @returns whether the link points to an existing volume accelerator
  template <concepts::accelerator_collection accel_coll_t,
            typename accel_index_t>
  DETRAY_HOST bool operator()(const accel_coll_t &coll,
                              const accel_index_t index) const {
    using accel_type = typename accel_coll_t::value_type;

    return concepts::volume_accelerator<accel_type> &&
           static_cast<dindex>(index) < coll.size();
  }
};

/// A functor that checks that the daughter volumes are inside their root
/// volume: All vertices of the daughter portals must be inside of the root
/// volume or on its boundary.
///
/// @note This is a necessary condition only: Curved portal boundaries are
/// tested at a finite number of points.
struct containment_checker {
  /// @tparam detector_t the type of the detector
  ///
  /// @param daughter_idx index of the daughter volume in the accelerator
  /// @param det the detector
  /// @param root_idx index of the root volume
  template <concepts::detector detector_t>
  DETRAY_HOST void operator()(const dindex daughter_idx, const detector_t &det,
                              const dindex root_idx) const {
    using algebra_t = typename detector_t::algebra_type;
    using scalar_t = dscalar<algebra_t>;
    using vertexer_t = vertexer<dpoint2D<algebra_t>, dpoint3D<algebra_t>>;

    // Tolerance for daughter portals on the boundary of the root volume
    constexpr scalar_t tol{10.f * unit<scalar_t>::um};
    // Number of segments for the curved portal boundaries
    constexpr dindex n_seg{16u};

    const typename detector_t::geometry_context ctx{};
    const auto root = tracking_volume{det, root_idx};

    for (const auto &pt_desc : tracking_volume{det, daughter_idx}.portals()) {
      const auto pt = detray::geometry::surface{det, pt_desc};
      const auto &trf = pt.transform(ctx);

      for (const auto &mask_vertices :
           pt.template visit_mask<vertexer_t>(n_seg)) {
        for (const auto &loc_vertex : mask_vertices) {
          const auto vertex = trf.point_to_global(loc_vertex);
          if (root.is_inside(vertex, ctx, tol) ==
              volume_containment::e_outside) {
            std::stringstream err_stream{};
            err_stream << "Daughter volume " << daughter_idx
                       << " is not contained in its mother volume " << root_idx
                       << ": Portal " << pt_desc.index()
                       << " has a vertex outside of the mother volume at ("
                       << vertex[0] << ", " << vertex[1] << ", " << vertex[2]
                       << ")";
            DETRAY_FATAL_HOST(err_stream.str());
            throw std::invalid_argument(err_stream.str());
          }
        }
      }
    }
  }
};

/// A functor that gets the vertices of the portal masks of a volume in local
/// coordinates (only for shapes that can form a volume boundary)
struct portal_vertices_getter {
  /// @tparam mask_group_t the type of the mask collection
  /// @tparam mask_range_t the type of the mask index (range) of the portal
  ///
  /// @param mask_group the collection of masks of the same type
  /// @param mask_range the index (range) of the masks of the portal
  /// @param n_seg the number of segments for curved mask boundaries
  ///
  /// @returns the vertices of all masks of the portal (local frame)
  template <typename mask_group_t, typename mask_range_t>
  DETRAY_HOST auto operator()(const mask_group_t &mask_group,
                              const mask_range_t &mask_range,
                              const dindex n_seg) const {
    using mask_t = typename mask_group_t::value_type;
    using shape_t = typename mask_t::shape;

    dvector<typename mask_t::point3_type> vertices{};
    if constexpr (std::same_as<shape_t, rectangle2D> ||
                  std::same_as<shape_t, trapezoid2D> ||
                  std::same_as<shape_t, annulus2D> ||
                  std::same_as<shape_t, ring2D> ||
                  std::same_as<shape_t, cylinder2D> ||
                  std::same_as<shape_t, concentric_cylinder2D>) {
      for (const auto &mask :
           detray::ranges::subrange(mask_group, mask_range)) {
        for (const auto &v : mask.vertices(n_seg)) {
          vertices.push_back(v);
        }
      }
    }
    return vertices;
  }
};

/// @brief Checks that all portal vertices of a volume lie on its boundary
///
/// The inside check of the volume uses only some of the portals (e.g. only
/// the discs and mantles of a cylinder). Make sure that the portals form a
/// consistent boundary, e.g. that no portal covers only a part of a face.
///
/// @param det the detector
/// @param vol the volume to be checked
/// @param err_stream the stream for error messages
///
/// @throws std::invalid_argument if a portal vertex is not on the boundary
template <concepts::detector detector_t, typename volume_t>
inline void check_portal_vertices(const detector_t &det, const volume_t &vol,
                                  std::stringstream &err_stream) {
  using scalar_t = dscalar<typename detector_t::algebra_type>;

  // Tolerance for the vertex positions
  constexpr scalar_t tol{10.f * unit<scalar_t>::um};
  // Number of segments for curved mask boundaries
  constexpr dindex n_seg{16u};

  const typename detector_t::geometry_context ctx{};

  for (const auto &pt_desc : vol.portals()) {
    const auto pt = detray::geometry::surface{det, pt_desc};
    const auto &trf = pt.transform(ctx);

    for (const auto &loc_vertex :
         pt.template visit_mask<portal_vertices_getter>(n_seg)) {
      const auto vertex = trf.point_to_global(loc_vertex);
      if (const auto vc = vol.is_inside(vertex, ctx, tol);
          vc != volume_containment::e_on_boundary) {
        err_stream << "Portal " << pt_desc.index() << " has a vertex that is "
                   << vc << " of its volume at (" << vertex[0] << ", "
                   << vertex[1] << ", " << vertex[2]
                   << "): The portals do not form a consistent volume "
                      "boundary:\n"
                   << vol;
        DETRAY_FATAL_HOST(err_stream.str());
        throw std::invalid_argument(err_stream.str());
      }
    }
  }
}

/// A functor that checks the material parametrization for a surface/volume
struct material_checker {
  /// Error message for material consistency check
  /// @param type Name of the material collection or object type.
  /// @param idx Index of the invalid material.
  /// @param mat Material value that failed validation.
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
    } catch (std::out_of_range & /*unused*/) {
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
                                     const index_t idx,
                                     const id_t /*unused*/) const {
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
        if (!mat.is_valid()) {
          throw_material_error("homogeneous surface material", idx, mat);
        }
      }
    } catch (std::out_of_range & /*unused*/) {
      std::stringstream err_stream{};
      err_stream << "Out of range material access in: "
                 << "homogeneous material collection at index " << idx;
      DETRAY_FATAL_HOST(err_stream.str());
      throw std::invalid_argument(err_stream.str());
    }
  }
};

/// @brief Checks the links between a volume and its mother and daughters
///
/// - The volume hierarchy has no cycles
/// - The accelerator link to the daughter volumes is valid
/// - The mother volume lists the volume exactly once as a daughter
/// - All daughters of the volume have the volume as mother
/// - All daughters are inside the volume (axis aligned bounding boxes)
///
/// @param det the detector
/// @param vol the volume to be checked
/// @param err_stream the stream for error messages
///
/// @throws std::invalid_argument if a check fails
template <concepts::detector detector_t, typename volume_t>
inline void check_volume_hierarchy(const detector_t &det, const volume_t &vol,
                                   std::stringstream &err_stream) {
  const auto n_volumes{static_cast<dindex>(det.volumes().size())};

  auto throw_error = [&err_stream]() {
    DETRAY_FATAL_HOST(err_stream.str());
    throw std::invalid_argument(err_stream.str());
  };

  // The daughter accelerator of a root volume must exist
  auto check_accel_link = [&det, &err_stream, &throw_error](const auto &v) {
    constexpr auto vol_obj_id{detector_t::geo_obj_ids::e_volume};
    if (!det.accelerator_store().template visit<volume_accel_link_checker>(
            det.volume(v.index()).template accel_link<vol_obj_id>())) {
      err_stream << "Invalid link to the daughter volume accelerator:\n" << v;
      throw_error();
    }
  };

  if (!vol.is_top_level()) {
    // Walk up the tree: Must reach a top-level volume
    dindex current{vol.mother()};
    dindex n_steps{0u};
    while (current != dindex_invalid) {
      if (current >= n_volumes) {
        err_stream << "Mother volume index " << current
                   << " out of bounds in the volume hierarchy of:\n"
                   << vol;
        throw_error();
      }
      if (current == vol.index() || n_steps >= n_volumes) {
        err_stream << "Cycle in volume hierarchy:\n" << vol;
        throw_error();
      }
      current = det.volume(current).mother();
      ++n_steps;
    }

    // The mother must list the volume exactly once
    const auto mother = tracking_volume{det, vol.mother()};
    dindex n_found{0u};
    if (mother.has_daughters()) {
      check_accel_link(mother);
      mother.template visit_daughter_volumes<daughter_checker>(
          det, mother.index(), vol.index(), n_found);
    }
    // TODO: Make this an error, once the volume grid builders are added
    if (n_found == 0u) {
      DETRAY_WARN_HOST(
          "Volume is not part of any volume accelerator of its "
          "mother volume "
          << mother.index() << ":\n"
          << vol);
    } else if (n_found != 1u) {
      err_stream << "Volume is registered " << n_found
                 << " times with its mother volume " << mother.index()
                 << " (expected once):\n"
                 << vol;
      throw_error();
    }
  }

  // All daughters must have this volume as mother
  if (vol.has_daughters()) {
    check_accel_link(vol);

    dindex n_found{0u};
    vol.template visit_daughter_volumes<daughter_checker>(
        det, vol.index(), dindex_invalid, n_found);

    vol.template visit_daughter_volumes<containment_checker>(det, vol.index());
  }
}

/// @brief Checks whether the data containers of a detector are empty
///
/// In case the default metadata is used, the unused containers are allowed to
/// be empty.
template <concepts::detector detector_t>
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
template <concepts::detector detector_t>
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

    // Check that the portals form a consistent volume boundary
    detail::check_portal_vertices(det, vol, err_stream);

    // Check the volume hierarchy
    detail::check_volume_hierarchy(det, vol, err_stream);

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

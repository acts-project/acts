// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/identifier.hpp"

// Sysytem include(s)
#include <memory>

namespace detray {

/// Templated surface class for detector surfaces and portals.
///
/// @note might be holding multiple surfaces in the future
///
/// @tparam mask_regsitry_t the type collection of masks that can be linked
///                         to the surface
/// @tparam material_registry_t the type collection of material that can be
///                             linked to the surface
/// @tparam transform_link_t how to reference the surfaces transforms
template <typename mask_link_t = dtyped_index<dindex, dindex>,
          typename material_link_t = dtyped_index<dindex, dindex>,
          typename transform_link_t = dindex,
          typename navigation_link_t = std::uint_least16_t>
class surface_descriptor {
 public:
  /// Link type of the mask to a volume.
  using navigation_link = navigation_link_t;
  // Broadcast the type of links
  using transform_link = transform_link_t;
  /// might be a single mask, a range of masks or a multiindex in the future
  using mask_link = mask_link_t;
  using mask_id = typename mask_link::id_type;
  using material_link = material_link_t;
  using material_id = typename material_link::id_type;

  /// Default constructor
  constexpr surface_descriptor() = default;

  /// Constructor with full arguments
  ///
  /// @param trf the transform for positioning and 3D local frame
  /// @param mask the type and index of the mask for this surface
  /// @param material the type and index of the material for this surface
  /// @param vol the volume this surface belongs to
  /// @param src the source object/source link this surface is representing
  /// @param sf_id remember whether this is a portal or not
  DETRAY_HOST
  constexpr surface_descriptor(const transform_link trf, const mask_link mask,
                               const material_link material,
                               const dindex volume, const surface_id sf_id)
      : m_identifier{geometry::identifier{}
                         .set_volume(volume)
                         .set_id(sf_id)
                         .set_transform(trf)},
        m_mask(mask),
        m_material(material) {}

  /// Equality operator
  ///
  /// @param rhs is the right hand side to be compared to
  DETRAY_HOST_DEVICE
  constexpr auto operator==(const surface_descriptor &rhs) const -> bool {
    return (m_mask == rhs.m_mask && m_material == rhs.m_material &&
            m_identifier == rhs.m_identifier);
  }

  /// Sets a new surface identifier
  DETRAY_HOST_DEVICE
  auto set_identifier(const geometry::identifier geo_id) -> void {
    m_identifier = geo_id;
  }

  /// @returns the surface identifier
  DETRAY_HOST_DEVICE
  constexpr auto identifier() const -> geometry::identifier {
    return m_identifier;
  }

  /// Sets a new surface id (portal/passive/sensitive)
  DETRAY_HOST_DEVICE
  auto set_id(const surface_id new_id) -> void { m_identifier.set_id(new_id); }

  /// @returns the surface id (sensitive, passive or portal)
  DETRAY_HOST_DEVICE
  constexpr auto id() const -> surface_id { return m_identifier.id(); }

  /// Sets a new volume link (index in volume collection of detector)
  DETRAY_HOST
  auto set_volume(const dindex new_idx) -> void {
    m_identifier.set_volume(new_idx);
  }

  /// @returns the surface id (sensitive, passive or portal)
  DETRAY_HOST_DEVICE
  constexpr auto volume() const -> dindex { return m_identifier.volume(); }

  /// Sets a new surface index (index in surface collection of surface store)
  DETRAY_HOST_DEVICE
  auto set_index(const dindex new_idx) -> void {
    m_identifier.set_index(new_idx);
  }

  /// @returns the surface id (sensitive, passive or portal)
  DETRAY_HOST_DEVICE
  constexpr auto index() const -> dindex { return m_identifier.index(); }

  /// Update the transform index
  ///
  /// @param offset update the position when move into new collection
  DETRAY_HOST
  auto update_transform(dindex offset) -> void {
    m_identifier.set_transform(transform() + offset);
  }

  /// @return the transform index
  DETRAY_HOST_DEVICE
  constexpr auto transform() const -> dindex {
    return m_identifier.transform();
  }

  /// Update the mask link
  ///
  /// @param offset update the position when move into new collection
  DETRAY_HOST
  auto update_mask(dindex offset) -> void { m_mask.shift(offset); }

  /// @return the mask link
  DETRAY_HOST_DEVICE
  constexpr auto mask() const -> const mask_link & { return m_mask; }

  /// Update the material link
  ///
  /// @param offset update the position when move into new collection
  DETRAY_HOST
  auto update_material(dindex offset) -> void { m_material.shift(offset); }

  /// Access to the material
  DETRAY_HOST_DEVICE
  constexpr auto material() -> material_link & { return m_material; }

  /// @return the material link
  DETRAY_HOST_DEVICE
  constexpr auto material() const -> const material_link & {
    return m_material;
  }

  /// @returns true if the surface descriptor has a valid material link
  DETRAY_HOST_DEVICE
  constexpr auto has_material() const -> bool {
    return (m_material.id() != material_link::id_type::e_none) &&
           !m_material.is_invalid();
  }

  /// @returns true if the surface is a sensitive detector module.
  DETRAY_HOST_DEVICE
  constexpr auto is_sensitive() const -> bool {
    return m_identifier.id() == surface_id::e_sensitive;
  }

  /// @returns true if the surface is a portal.
  DETRAY_HOST_DEVICE
  constexpr auto is_portal() const -> bool {
    return m_identifier.id() == surface_id::e_portal;
  }

  /// @returns true if the surface is a passive detector element.
  DETRAY_HOST_DEVICE
  constexpr auto is_passive() const -> bool {
    return m_identifier.id() == surface_id::e_passive;
  }

  /// @returns a string stream that prints the surface details
  DETRAY_HOST
  friend std::ostream &operator<<(std::ostream &os,
                                  const surface_descriptor &sf) {
    os << sf.m_identifier;
    os << " | mask: " << sf.m_mask;
    os << " | mat.: " << sf.m_material;
    return os;
  }

 private:
  geometry::identifier m_identifier{};
  mask_link m_mask{};
  material_link m_material{};
};

}  // namespace detray

namespace std {

// Specialize std::hash so surface descriptors can be used in STL containers
template <typename mask_link_t, typename material_link_t,
          typename transform_link_t, typename navigation_link_t>
struct hash<detray::surface_descriptor<mask_link_t, material_link_t,
                                       transform_link_t, navigation_link_t>> {
  using descr_t =
      detray::surface_descriptor<mask_link_t, material_link_t, transform_link_t,
                                 navigation_link_t>;

  auto operator()(const descr_t sf_desc) const noexcept {
    return std::hash<detray::geometry::identifier::value_t>()(
        sf_desc.identifier().value());
  }
};

}  // namespace std

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/identifier.hpp"
#include "detray/geometry/mask.hpp"

namespace detray {

template <typename algebra_t, typename shape_t, typename material_t,
          typename nav_link_t>
class external_surface {
 public:
  using transform3_type = dtransform3D<algebra_t>;
  using mask_type = detray::mask<shape_t, algebra_t, nav_link_t>;
  using material_type = material_t;
  /// Link type of the mask to a volume.
  using navigation_link = nav_link_t;

  /// Default constructor
  constexpr external_surface() = default;

  /// Constructor with full arguments
  ///
  /// @param trf the transform for positioning and 3D local frame
  /// @param mask the surface mask
  /// @param material the surface material
  /// @param volume the volume this surface belongs to
  /// @param sf_id remember whether this is a portal or not
  DETRAY_HOST_DEVICE
  constexpr external_surface(const transform3_type& trf, const mask_type& mask,
                             const material_type& material, const dindex volume,
                             const surface_id sf_id)
      : m_transform{trf},
        m_mask(mask),
        m_material(material),
        m_barcode{geometry::identifier{}.set_volume(volume).set_id(sf_id)} {}

  /// Constructor with full arguments - move semantics
  ///
  /// @param trf the transform for positioning and 3D local frame
  /// @param mask the surface mask
  /// @param material the surface material
  /// @param volume the volume this surface belongs to
  /// @param sf_id remember whether this is a portal or not
  DETRAY_HOST_DEVICE
  constexpr external_surface(transform3_type&& trf, mask_type&& mask,
                             material_type&& material, const dindex volume,
                             const surface_id sf_id)
      : m_transform{std::move(trf)},
        m_mask(std::move(mask)),
        m_material(std::move(material)),
        m_barcode{geometry::identifier{}.set_volume(volume).set_id(sf_id)} {}

  /// Equality operator
  ///
  /// @param rhs is the right hand side to be compared to
  DETRAY_HOST_DEVICE
  constexpr auto operator==(const external_surface& rhs) const -> bool {
    return (m_barcode == rhs.m_barcode);
  }

  /// Sets a new surface barcode
  DETRAY_HOST_DEVICE
  auto set_barcode(const geometry::identifier bcd) -> void { m_barcode = bcd; }

  /// @returns the surface barcode
  DETRAY_HOST_DEVICE
  constexpr auto barcode() const -> geometry::identifier { return m_barcode; }

  /// @returns the surface identifier (identifier-era naming)
  DETRAY_HOST_DEVICE
  constexpr auto identifier() const -> geometry::identifier {
    return m_barcode;
  }

  /// Explicit conversion for APIs that take geometry::identifier.
  DETRAY_HOST_DEVICE
  explicit operator geometry::identifier() const {
    return m_barcode;
  }

  /// Sets a new surface id (portal/passive/sensitive)
  DETRAY_HOST_DEVICE
  auto set_id(const surface_id new_id) -> void { m_barcode.set_id(new_id); }

  /// @returns the surface id (sensitive, passive or portal)
  DETRAY_HOST_DEVICE
  constexpr auto id() const -> surface_id { return m_barcode.id(); }

  /// Sets a new volume link (index in volume collection of detector)
  DETRAY_HOST_DEVICE
  auto set_volume(const dindex new_idx) -> void {
    m_barcode.set_volume(new_idx);
  }

  /// @returns the surface id (sensitive, passive or portal)
  DETRAY_HOST_DEVICE
  constexpr auto volume() const -> dindex { return m_barcode.volume(); }

  /// Sets a new surface index (index in surface collection of surface store)
  DETRAY_HOST_DEVICE
  auto set_index(const dindex new_idx) -> void { m_barcode.set_index(new_idx); }

  /// @returns the surface id (sensitive, passive or portal)
  DETRAY_HOST_DEVICE
  constexpr auto index() const -> dindex { return m_barcode.index(); }

  /// Update the transform index
  ///
  /// @param offset update the position when move into new collection
  DETRAY_HOST_DEVICE
  void set_transform(const transform3_type& trf) { m_transform = trf; }

  /// @return the transform index
  DETRAY_HOST_DEVICE
  constexpr auto transform() const -> const transform3_type& {
    return m_transform;
  }

  /// Update the mask link
  ///
  /// @param offset update the position when move into new collection
  DETRAY_HOST_DEVICE
  void set_mask(const mask_type m) { m_mask = m; }

  /// @return the mask - const
  DETRAY_HOST_DEVICE
  constexpr auto mask() const -> const mask_type& { return m_mask; }

  /// Update the material link
  ///
  /// @param offset update the position when move into new collection
  DETRAY_HOST_DEVICE
  void set_material(const material_type& offset) { m_material.shift(offset); }

  /// Access to the material
  DETRAY_HOST_DEVICE
  constexpr auto material() -> material_type& { return m_material; }

  /// @return the material link
  DETRAY_HOST_DEVICE
  constexpr auto material() const -> const material_type& { return m_material; }

  /// @returns true if the surface descriptor has a valid material link
  DETRAY_HOST_DEVICE
  constexpr auto has_material() const -> bool { return false; }

  /// @returns true if the surface is a sensitive detector module.
  DETRAY_HOST_DEVICE
  constexpr auto is_sensitive() const -> bool {
    return m_barcode.id() == surface_id::e_sensitive;
  }

  /// @returns true if the surface is a portal.
  DETRAY_HOST_DEVICE
  constexpr auto is_portal() const -> bool {
    return m_barcode.id() == surface_id::e_portal;
  }

  /// @returns true if the surface is a passive detector element.
  DETRAY_HOST_DEVICE
  constexpr auto is_passive() const -> bool {
    return m_barcode.id() == surface_id::e_passive;
  }

 private:
  /// Surface placement transform and its inverse
  transform3_type m_transform{};
  /// Surface mask
  mask_type m_mask{};
  /// Surface material
  material_type m_material{};
  /// Surface hash
  geometry::identifier m_barcode{};
};

}  // namespace detray

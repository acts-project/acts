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
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/math.hpp"
#include "detray/geometry/detail/surface_kernels.hpp"
#include "detray/geometry/identifier.hpp"
#include "detray/material/material.hpp"
#include "detray/utils/ranges/ranges.hpp"

// System include(s)
#include <ostream>
#include <type_traits>

namespace detray::geometry {

/// @brief Facade for a detray detector surface.
///
/// Provides an interface to geometry specific functionality like
/// local-to-global coordinate transforms or mask and material visitors. It
/// wraps a detector instance that contains the data and a surface descriptor
/// that contains the indices into the detector data containers for the
/// specific surface instance.
template <typename det_t>  // @TODO: This needs a concept
class surface {
  /// Make sure the detector is always evaluated as constant type
  using detector_t = std::add_const_t<det_t>;
  /// Surface descriptor type
  using descr_t = typename detector_t::surface_type;
  /// Implementation
  using kernels = detail::surface_kernels<typename detector_t::algebra_type>;

 public:
  using algebra_type = typename detector_t::algebra_type;
  using scalar_type = dscalar<algebra_type>;
  using point2_type = dpoint2D<algebra_type>;
  using point3_type = dpoint3D<algebra_type>;
  using vector3_type = dvector3D<algebra_type>;
  using transform3_type = dtransform3D<algebra_type>;
  using context = typename detector_t::geometry_context;

  /// Not allowed: always needs a detector and a descriptor.
  surface() = delete;

  /// Constructor from detector @param det and surface descriptor
  /// @param desc from that detector.
  DETRAY_HOST_DEVICE
  constexpr surface(const detector_t &det, const descr_t &desc)
      : m_detector{det}, m_desc{desc} {
    assert(!m_desc.identifier().is_invalid());
    assert(m_desc.index() < det.surfaces().size());
    assert(m_desc.transform() < det.transform_store().size());
  }

  /// Constructor from detector @param det and identifier @param geo_id in
  /// that detector.
  DETRAY_HOST_DEVICE
  constexpr surface(const detector_t &det, const geometry::identifier geo_id)
      : surface(det, det.surface(geo_id)) {}

  /// Constructor from detector @param det and surface index @param sf_idx
  DETRAY_HOST_DEVICE
  constexpr surface(const detector_t &det, const dindex sf_idx)
      : surface(det, det.surface(sf_idx)) {}

  /// Allow conversion from surface<detector_t> to surface<const detector_t>,
  /// which has no semantic difference (the detector is always const internally)
  template <typename detector_type = detector_t>
    requires(!std::is_const_v<detector_type>)
  // NOLINTNEXTLINE
  DETRAY_HOST_DEVICE constexpr operator surface<const detector_type>() const {
    return surface<const detector_type>{this->m_detector, this->m_desc};
  }

  /// Equality operator
  ///
  /// @param rhs is the right hand side to be compared to
  DETRAY_HOST_DEVICE
  constexpr auto operator==(const surface &rhs) const -> bool {
    return (&m_detector == &(rhs.m_detector) && m_desc == rhs.m_desc);
  }

  /// @returns the surface identifier
  DETRAY_HOST_DEVICE
  constexpr auto identifier() const -> geometry::identifier {
    assert(!m_desc.identifier().is_invalid());
    return m_desc.identifier();
  }

  /// @returns the index of the mother volume
  DETRAY_HOST_DEVICE
  constexpr auto volume() const -> dindex {
    assert(identifier().volume() < m_detector.volumes().size());
    return identifier().volume();
  }

  /// @returns the index of the surface in the detector surface lookup
  DETRAY_HOST_DEVICE
  constexpr auto index() const -> dindex {
    assert(identifier().index() < m_detector.surfaces().size());
    return identifier().index();
  }

  /// @returns the surface id (sensitive, passive or portal)
  DETRAY_HOST_DEVICE
  constexpr auto id() const -> surface_id {
    assert(identifier().id() != surface_id::e_unknown);
    return identifier().id();
  }

  /// @returns the extra bits in the identifier
  DETRAY_HOST_DEVICE
  constexpr auto extra() const -> dindex { return identifier().extra(); }

  /// @returns an id for the surface type (e.g. 'rectangle')
  DETRAY_HOST_DEVICE
  constexpr auto shape_id() const { return m_desc.mask().id(); }

  /// @returns the surface source link
  DETRAY_HOST_DEVICE
  constexpr auto source() const {
    return m_detector.surface(identifier()).source;
  }

  /// @returns true if the surface is a sensitive detector module.
  DETRAY_HOST_DEVICE
  constexpr auto is_sensitive() const -> bool {
    return identifier().id() == surface_id::e_sensitive;
  }

  /// @returns true if the surface is a portal.
  DETRAY_HOST_DEVICE
  constexpr auto is_portal() const -> bool {
    return identifier().id() == surface_id::e_portal;
  }

  /// @returns true if the surface is a passive detector element.
  DETRAY_HOST_DEVICE
  constexpr auto is_passive() const -> bool {
    return identifier().id() == surface_id::e_passive;
  }

  /// @returns true if the surface carries material.
  DETRAY_HOST_DEVICE
  constexpr bool has_material() const { return m_desc.has_material(); }

  /// @returns the mask volume link
  DETRAY_HOST_DEVICE
  constexpr auto volume_links() const {
    return visit_mask<typename kernels::get_volume_links>();
  }

  /// @returns the mask shape name
  DETRAY_HOST
  std::string shape_name() const {
    return visit_mask<typename kernels::get_shape_name>();
  }

  /// @returns the number of masks associated with this surface
  DETRAY_HOST_DEVICE
  std::size_t n_masks() const {
    if constexpr (concepts::interval<decltype(m_desc.mask().index())>) {
      return m_desc.mask().index().size();
    } else {
      return 1u;
    }
  }

  /// @returns the coordinate transform matrix of the surface
  DETRAY_HOST_DEVICE
  constexpr auto transform(const context &ctx) const
      -> const transform3_type & {
    assert(m_desc.transform() < m_detector.transform_store().size());
    return m_detector.transform_store().at(m_desc.transform(), ctx);
  }

  /// @returns the mask volume link
  template <concepts::point point_t = point2_type>
  DETRAY_HOST_DEVICE constexpr bool is_inside(const point_t &loc_p,
                                              const scalar_type tol) const {
    return visit_mask<typename kernels::is_inside>(loc_p, tol);
  }

  /// @returns a boundary value of the surface, according to @param index
  DETRAY_HOST_DEVICE
  constexpr auto boundary(std::size_t index) const {
    return visit_mask<typename kernels::get_mask_value>(index);
  }

  /// @returns the centroid of the surface mask in local cartesian coordinates.
  DETRAY_HOST_DEVICE
  constexpr auto centroid() const -> point3_type {
    return visit_mask<typename kernels::centroid>();
  }

  /// @returns the center position of the surface in global coordinates
  /// @note for shapes like the annulus this is not synonymous to the centroid
  /// but instead the focal point of the strip system
  DETRAY_HOST_DEVICE
  constexpr auto center(const context &ctx) const -> point3_type {
    return transform(ctx).translation();
  }

  /// @returns the surface normal in global coordinates at a given bound/local
  /// position @param p
  template <concepts::point point_t = point2_type>
  DETRAY_HOST_DEVICE constexpr auto normal(const context &ctx,
                                           const point_t &p) const
      -> vector3_type {
    return visit_mask<typename kernels::normal>(transform(ctx), p);
  }

  /// @returns a pointer to the material parameters at the local position
  /// @param loc_p
  DETRAY_HOST_DEVICE constexpr const material<scalar_type> *material_parameters(
      const point2_type &loc_p) const {
    return visit_material<typename kernels::get_material_params>(loc_p);
  }

  /// @returns the local position to the global point @param global for
  /// a given geometry context @param ctx and track direction @param dir
  DETRAY_HOST_DEVICE
  constexpr point3_type global_to_local(const context &ctx,
                                        const point3_type &global,
                                        const vector3_type &dir) const {
    return visit_mask<typename kernels::global_to_local>(transform(ctx), global,
                                                         dir);
  }

  /// @returns the bound (2D) position to the global point @param global for
  /// a given geometry context @param ctx and track direction @param dir
  DETRAY_HOST_DEVICE
  constexpr point2_type global_to_bound(const context &ctx,
                                        const point3_type &global,
                                        const vector3_type &dir) const {
    const point3_type local{global_to_local(ctx, global, dir)};

    return {local[0], local[1]};
  }

  /// @returns the global position to the given local position @param local
  /// for a given geometry context @param ctx
  template <concepts::point point_t = point2_type>
  DETRAY_HOST_DEVICE constexpr point3_type local_to_global(
      const context &ctx, const point_t &local, const vector3_type &dir) const {
    return visit_mask<typename kernels::local_to_global>(transform(ctx), local,
                                                         dir);
  }

  /// Call a functor on the surfaces mask with additional arguments.
  ///
  /// @tparam functor_t the prescription to be applied to the mask
  /// @tparam Args      types of additional arguments to the functor
  template <typename functor_t, typename... Args>
  DETRAY_HOST_DEVICE constexpr auto visit_mask(Args &&...args) const {
    assert(!m_desc.mask().is_invalid());
    const auto &masks = m_detector.mask_store();
    return masks.template visit<functor_t>(m_desc.mask(),
                                           std::forward<Args>(args)...);
  }

  /// Call a functor on the surfaces material with additional arguments.
  ///
  /// @tparam functor_t the prescription to be applied to the material
  /// @tparam Args      types of additional arguments to the functor
  template <typename functor_t, typename... Args>
  DETRAY_HOST_DEVICE constexpr auto visit_material(Args &&...args) const {
    assert(has_material());
    const auto &materials = m_detector.material_store();
    return materials.template visit<functor_t>(m_desc.material(),
                                               std::forward<Args>(args)...);
  }

  /// Do a consistency check on the surface after building the detector.
  ///
  /// @param os output stream for error messages.
  ///
  /// @returns true if the surface is consistent
  DETRAY_HOST bool self_check(std::ostream &os) const {
    if (identifier().is_invalid()) {
      os << "DETRAY ERROR (HOST): Invalid identifier for surface:\n"
         << *this << std::endl;
      return false;
    }
    if (index() >= m_detector.surfaces().size()) {
      os << "DETRAY ERROR (HOST): Surface index out of bounds for "
            "surface:\n"
         << *this << std::endl;
      return false;
    }
    if (volume() >= m_detector.volumes().size()) {
      os << "DETRAY ERROR (HOST): Surface volume index out of bounds for "
            "surface:\n"
         << *this << std::endl;
      return false;
    }
    if (detail::is_invalid_value(m_desc.transform())) {
      os << "DETRAY ERROR (HOST): Surface transform undefined for "
            "surface:\n"
         << *this << std::endl;
      return false;
    }
    if (m_desc.transform() >= m_detector.transform_store().size()) {
      os << "DETRAY ERROR (HOST): Surface transform index out of bounds "
            "for surface:\n"
         << *this << std::endl;
      return false;
    }
    if (detail::is_invalid_value(m_desc.mask())) {
      os << "DETRAY ERROR (HOST): Surface does not have a valid mask "
            "link:\n"
         << *this << std::endl;
      return false;
    }
    // Only check, if there is material in the detector
    if (!m_detector.material_store().all_empty() && has_material() &&
        m_desc.material().is_invalid_index()) {
      os << "DETRAY ERROR (HOST): Surface does not have valid material "
            "link:\n"
         << *this << std::endl;
      return false;
    }
    // Check the mask boundaries
    if (!visit_mask<typename kernels::mask_self_check>(os)) {
      os << "\nSurface: " << *this << std::endl;
      return false;
    }

    // Check the mask volume link
    const auto vol_links = visit_mask<typename kernels::get_volume_links>();
    for (const auto vol_link : vol_links) {
      if (is_portal()) {
        if (vol_link == volume()) {
          os << "DETRAY ERROR (HOST): Portal surface links to mother "
                "volume:\n"
             << *this << std::endl;
          return false;
        }
      } else if (vol_link != volume()) {
        os << "DETRAY ERROR (HOST): Passive/sensitive surface does not "
              "link to mother volume: Mask volume link : "
           << vol_link << "\n"
           << *this << std::endl;
        return false;
      }
    }

    return true;
  }

 protected:
  /// @returns a string stream that prints the surface details
  DETRAY_HOST
  friend std::ostream &operator<<(std::ostream &os, const surface &sf) {
    os << sf.m_desc;
    return os;
  }

  /// @returns access to the underlying detector object
  DETRAY_HOST_DEVICE
  const detector_t &detector() const { return m_detector; }

  /// @returns access to the underlying surface descriptor object
  DETRAY_HOST_DEVICE
  descr_t descriptor() const { return m_desc; }

 private:
  /// Access to the detector stores
  const detector_t &m_detector;
  /// Access to the descriptor
  const descr_t m_desc;
};

template <typename detector_t, typename descr_t>
DETRAY_HOST_DEVICE surface(const detector_t &, const descr_t &)
    -> surface<detector_t>;

template <typename detector_t>
DETRAY_HOST_DEVICE surface(const detector_t &, const geometry::identifier)
    -> surface<detector_t>;

}  // namespace detray::geometry

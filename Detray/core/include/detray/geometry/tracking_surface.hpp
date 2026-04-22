// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/geometry/detail/tracking_surface_kernels.hpp"
#include "detray/geometry/surface.hpp"

namespace detray {

/// @brief Facade for a detray detector surface with extra tracking capabilities
template <typename detector_t>  // @TODO: This needs a concept
class tracking_surface : public geometry::surface<detector_t> {
  using base_surface_t = geometry::surface<detector_t>;

  /// Implementation of tracking functionality
  using kernels =
      detail::tracking_surface_kernels<typename detector_t::algebra_type>;
  /// Vector type for track parameters in global coordinates
  using free_param_vector_type = typename kernels::free_param_vector_type;
  /// Vector type for track parameters in local (bound) coordinates
  using bound_param_vector_type = typename kernels::bound_param_vector_type;

 public:
  using algebra_type = typename base_surface_t::algebra_type;
  using scalar_type = typename base_surface_t::scalar_type;
  using point2_type = typename base_surface_t::point2_type;
  using point3_type = typename base_surface_t::point3_type;
  using vector3_type = typename base_surface_t::vector3_type;
  using transform3_type = typename base_surface_t::transform3_type;
  using context = typename base_surface_t::context;

  /// Use base class constructors
  using base_surface_t::base_surface_t;

  /// Decorate a geometric surface with tracking functionality
  DETRAY_HOST_DEVICE
  explicit constexpr tracking_surface(const base_surface_t sf)
      : base_surface_t(sf) {}

  /// Conversion to surface interface around constant detector type
  template <typename detector_type = detector_t>
    requires(!std::is_const_v<detector_type>)
  DETRAY_HOST_DEVICE constexpr operator tracking_surface<const detector_type>()
      const {
    return tracking_surface<const detector_type>{this->m_detector,
                                                 this->m_desc};
  }

  /// @returns the track parametrization projected onto the surface (bound)
  DETRAY_HOST_DEVICE
  constexpr auto free_to_bound_vector(
      const context &ctx, const free_param_vector_type &free_vec) const {
    return this->template visit_mask<typename kernels::free_to_bound_vector>(
        this->transform(ctx), free_vec);
  }

  /// @returns the global track parametrization from a bound representation
  DETRAY_HOST_DEVICE
  constexpr auto bound_to_free_vector(
      const context &ctx, const bound_param_vector_type &bound_vec) const {
    return this->template visit_mask<typename kernels::bound_to_free_vector>(
        this->transform(ctx), bound_vec);
  }

  /// @returns the jacobian to go from a free to a bound track parametrization
  DETRAY_HOST_DEVICE
  constexpr auto free_to_bound_jacobian(
      const context &ctx, const free_param_vector_type &free_vec) const {
    return this->template visit_mask<typename kernels::free_to_bound_jacobian>(
        this->transform(ctx), free_vec);
  }

  /// @returns the jacobian to go from a bound to a free track parametrization
  DETRAY_HOST_DEVICE
  constexpr auto bound_to_free_jacobian(
      const context &ctx, const bound_param_vector_type &bound_vec) const {
    return this->template visit_mask<typename kernels::bound_to_free_jacobian>(
        this->transform(ctx), bound_vec);
  }

  /// @returns the path correction term
  DETRAY_HOST_DEVICE
  constexpr auto path_correction(const context &ctx, const vector3_type &pos,
                                 const vector3_type &dir,
                                 const vector3_type &dtds,
                                 const scalar_type dqopds) const {
    return this->template visit_mask<typename kernels::path_correction>(
        this->transform(ctx), pos, dir, dtds, dqopds);
  }

  /// @returns the free to path derivative
  DETRAY_HOST_DEVICE
  constexpr auto free_to_path_derivative(const context &ctx,
                                         const vector3_type &pos,
                                         const vector3_type &dir,
                                         const vector3_type &dtds) const {
    return this->template visit_mask<typename kernels::free_to_path_derivative>(
        this->transform(ctx), pos, dir, dtds);
  }
};

template <typename detector_t, typename descr_t>
DETRAY_HOST_DEVICE tracking_surface(const detector_t &, const descr_t &)
    -> tracking_surface<detector_t>;

template <typename detector_t>
DETRAY_HOST_DEVICE tracking_surface(const detector_t &,
                                    const geometry::identifier)
    -> tracking_surface<detector_t>;

template <typename detector_t>
DETRAY_HOST_DEVICE tracking_surface(const geometry::surface<detector_t>)
    -> tracking_surface<detector_t>;

}  // namespace detray

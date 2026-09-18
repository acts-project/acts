/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <detray/definitions/algebra.hpp>
#include <detray/geometry/identifier.hpp>
#include <traccc/definitions/qualifiers.hpp>
#include <traccc/edm/track_parameters.hpp>
#include <traccc/utils/soa.hpp>

namespace traccc {

template <typename algebra_t>
using soa_bound_track_parameters_layout =
    soa_layout<detray::dscalar<algebra_t>[6], detray::dscalar<algebra_t>[36],
               detray::geometry::identifier>;

template <typename algebra_t>
using soa_bound_track_parameters_buffer =
    soa_buffer<soa_bound_track_parameters_layout<algebra_t>>;

template <typename algebra_t>
using soa_bound_track_parameters_view =
    soa_view<soa_bound_track_parameters_layout<algebra_t>>;

template <typename algebra_t>
using soa_bound_track_parameters_const_view =
    soa_const_view<soa_bound_track_parameters_layout<algebra_t>>;

template <typename algebra_t, bool is_const = false>
struct soa_bound_track_parameters_device
    : private soa_device<soa_bound_track_parameters_layout<algebra_t>, true,
                         is_const> {
  using base_type =
      soa_device<soa_bound_track_parameters_layout<algebra_t>, true, is_const>;
  using base_type::base_type;
  using value_type = bound_track_parameters<algebra_t>;
  using scalar_type = detray::dscalar<algebra_t>;

  TRACCC_HOST_DEVICE value_type get(unsigned int x) const {
    value_type v;
    typename value_type::covariance_type cov;

    for (unsigned int i = 0; i < 6; ++i) {
      v[i] = vector(x, i);
    }

    for (unsigned int i = 0; i < 6; ++i) {
      for (unsigned int j = 0; j < 6; ++j) {
        detray::getter::element(cov, i, j) = covariance(x, i * 6 + j);
      }
    }

    v.set_covariance(cov);
    v.set_surface_link(surface_link(x));

    return v;
  }

  TRACCC_HOST_DEVICE void set(unsigned int x, const value_type& v)
    requires(!is_const)
  {
    surface_link(x) = v.surface_link();

    for (unsigned int i = 0; i < 6; ++i) {
      vector(x, i) = v[i];
    }

    for (unsigned int i = 0; i < 6; ++i) {
      for (unsigned int j = 0; j < 6; ++j) {
        covariance(x, i * 6 + j) =
            detray::getter::element(v.covariance(), i, j);
      }
    }
  }

  TRACCC_HOST_DEVICE scalar_type& vector(unsigned int x, unsigned int i)
    requires(!is_const)
  {
    return this->template at<0>(x, i);
  }

  TRACCC_HOST_DEVICE scalar_type vector(unsigned int x, unsigned int i) const {
    return this->template at<0>(x, i);
  }

  TRACCC_HOST_DEVICE scalar_type& covariance(unsigned int x, unsigned int i)
    requires(!is_const)
  {
    return this->template at<1>(x, i);
  }

  TRACCC_HOST_DEVICE scalar_type covariance(unsigned int x,
                                            unsigned int i) const {
    return this->template at<1>(x, i);
  }

  TRACCC_HOST_DEVICE detray::geometry::identifier& surface_link(unsigned int x)
    requires(!is_const)
  {
    return this->template at<2>(x);
  }

  TRACCC_HOST_DEVICE detray::geometry::identifier surface_link(
      unsigned int x) const {
    return this->template at<2>(x);
  }
};

template <typename algebra_t>
using soa_bound_track_parameters_const_device =
    soa_bound_track_parameters_device<algebra_t, true>;

}  // namespace traccc

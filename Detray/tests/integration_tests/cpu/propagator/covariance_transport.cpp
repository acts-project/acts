// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/definitions/units.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes.hpp"
#include "detray/geometry/shapes/unbounded.hpp"
#include "detray/geometry/surface_descriptor.hpp"
#include "detray/navigation/intersection/helix_intersector.hpp"
#include "detray/propagator/detail/jacobian_engine.hpp"
#include "detray/test/framework/types.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/axis_rotation.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// google-test include(s).
#include <gtest/gtest.h>

using namespace detray;

// Algebra types
using test_algebra = test::algebra;
using scalar = test::scalar;
using transform3 = test::transform3;
using vector3 = test::vector3;
using intersection_t = intersection2D<surface_descriptor<>, test_algebra,
                                      intersection::contains_pos>;

// Mask types to be tested
// @TODO: Remove unbounded tag
using annulus_type =
    detray::mask<detray::unbounded<detray::annulus2D>, test_algebra>;
using rectangle_type = detray::mask<detray::rectangle2D, test_algebra>;
using trapezoid_type = detray::mask<detray::trapezoid2D, test_algebra>;
using ring_type = detray::mask<detray::ring2D, test_algebra>;
using cylinder_type = detray::mask<detray::cylinder2D, test_algebra>;
using straw_tube_type = detray::mask<detray::line_circular, test_algebra>;
using drift_cell_type = detray::mask<detray::line_square, test_algebra>;

// Test class for covariance transport
template <typename T>
class detray_propagation_HelixCovarianceTransportValidation
    : public ::testing::Test {
 public:
  // Environment Setup
  const vector3 B{static_cast<scalar>(0.f), static_cast<scalar>(0.f),
                  1.f * unit<scalar>::T};
  const vector3 y_axis{0.f, 1.f, 0.f};
  const vector3 z_axis{0.f, 0.f, 1.f};
  static constexpr const scalar mask_tolerance{1e-3f};
  static constexpr const scalar tolerance{3e-2f};

  // Test types
  using mask_type = T;
  using local_frame_type = typename mask_type::local_frame;

  // First mask at the origin is always rectangle
  using first_mask_type = rectangle_type;
  using first_local_frame_type = typename first_mask_type::local_frame;

  // Algebra type
  using algebra_type = typename local_frame_type::algebra_type;

  // Vector and matrix types
  using bound_param_vector_t =
      typename bound_track_parameters<algebra_type>::parameter_vector_type;
  using bound_matrix_t = bound_matrix<algebra_type>;
  using bound_to_free_matrix_t = bound_to_free_matrix<algebra_type>;

  using free_matrix_t = free_matrix<algebra_type>;
  using free_to_bound_matrix_t = free_to_bound_matrix<algebra_type>;

  std::tuple<annulus_type, rectangle_type, trapezoid_type, ring_type,
             cylinder_type, straw_tube_type, drift_cell_type>
      masks = std::make_tuple(
          annulus_type{0u, 7.2f * unit<scalar>::mm, 12.0f * unit<scalar>::mm,
                       0.74195f, 1.33970f, 0.f, -2.f, 2.f},
          rectangle_type{0u, 50 * unit<scalar>::mm, 50 * unit<scalar>::mm},
          trapezoid_type{0u, 50 * unit<scalar>::mm, 100 * unit<scalar>::mm,
                         50 * unit<scalar>::mm,
                         1.f / (2.f * 50 * unit<scalar>::mm)},
          ring_type{0u, 0.f, 50.f * unit<scalar>::mm},
          cylinder_type{0u, 0.3f * unit<scalar>::mm, -100.f * unit<scalar>::mm,
                        100.f * unit<scalar>::mm},
          straw_tube_type{0u, 50.f * unit<scalar>::mm,
                          100.f * unit<scalar>::mm},
          drift_cell_type{0u, 50.f * unit<scalar>::mm,
                          100.f * unit<scalar>::mm});

  // Create transform matrices
  std::vector<transform3> create_transforms(
      detail::helix<algebra_type>& reference_helix,
      const std::size_t n_planes) {
    std::vector<transform3> trfs;

    // Step size between two neighbor planes
    const scalar S{2.f * constant<scalar>::pi /
                   std::abs(reference_helix.qop() * reference_helix.B())};
    const scalar step_size = S / static_cast<scalar>(n_planes);

    for (std::size_t i = 0u; i < n_planes; i++) {
      const scalar s = step_size * scalar(i);

      // Translation of the new surface
      vector3 trl = reference_helix(s);

      // Normal vector of the new surface
      vector3 w = reference_helix.dir(s);
      vector3 v = vector::cross(z_axis, w);

      if (i > 0u &&
          (std::is_same_v<local_frame_type, cylindrical2D<algebra_type>> ||
           std::is_same_v<local_frame_type, line2D<algebra_type>>)) {
        const vector3 r_axis = vector::cross(w, z_axis);

        // Rotate for cylinder and wire
        axis_rotation<algebra_type> axis_rot(r_axis,
                                             constant<scalar>::pi / 2.f);

        // Test masks are rotated
        w = axis_rot(w);
        EXPECT_NEAR(vector::norm(r_axis), 1.f, tolerance);

        v = vector::cross(r_axis, w);

        /*
        const vector3 offset{0.f, 0.f, 10.f * unit<scalar>::mm};
        trl = trl + offset;
        */
      } else {
        // @TODO why does this offset (in y-direction) fail the test???
        // const vector3 offset{0.f, 10.f * unit<scalar>::mm, 10.f *
        // unit<scalar>::mm};
        const vector3 offset{static_cast<scalar>(0.f), static_cast<scalar>(0.f),
                             10.f * unit<scalar>::mm};
        trl = trl + offset;
      }

      // Add transform matrix
      trfs.emplace_back(trl, w, v);
    }

    return trfs;
  }

  // Error propagation
  template <typename departure_mask_type, typename destination_mask_type>
  bound_track_parameters<algebra_type> propagate(
      const bound_track_parameters<algebra_type>& bound_params,
      const transform3& trf_0, const transform3& trf_1,
      const departure_mask_type& mask_0, const destination_mask_type& mask_1,
      scalar& total_path_length, std::vector<intersection_t>& sfis) {
    using departure_frame = typename departure_mask_type::local_frame;
    using destination_frame = typename destination_mask_type::local_frame;

    using jacobian_engine = detail::jacobian_engine<algebra_type>;

    const bound_param_vector_t& bound_vec_0 = bound_params;
    const bound_matrix_t& bound_cov_0 = bound_params.covariance();

    // Free vector at the departure surface
    const auto free_trk_0 =
        detail::bound_to_free_vector(trf_0, mask_0, bound_vec_0);

    // Helix from the departure surface
    detail::helix<algebra_type> hlx(free_trk_0, B);

    // Bound-to-free jacobian at the departure surface
    const bound_to_free_matrix_t bound_to_free_jacobi =
        jacobian_engine::template bound_to_free_jacobian<departure_frame>(
            trf_0, mask_0, bound_vec_0);

    // Get the intersection on the next surface
    helix_intersector<typename destination_mask_type::shape, algebra_type>
        helix_inters{};
    helix_inters.run_rtsafe = false;
    const intersection_t is = get_intersection(helix_inters(
        hlx, surface_descriptor<>{}, mask_1, trf_1, this->mask_tolerance));
    // Check for successful intersection
    EXPECT_TRUE(is.is_inside()) << is;

    sfis.push_back(is);

    // Helical path length between two surfaces
    const auto path_length = is.path();

    // Add the path length to the total path length
    total_path_length += path_length;

    // Free transport jacobian between two surfaces
    const free_matrix_t transport_jacobi = hlx.jacobian(path_length);

    // r at the destination surface
    const vector3 r = hlx.pos(path_length);

    // dr/ds, or t at the destination surface
    const vector3 t = hlx.dir(path_length);

    // d^2r/ds^2, or dt/ds at the destination surface
    const vector3 dtds = hlx.qop() * vector::cross(t, B);

    // d(qop)/ds, which is zero in this test
    const scalar dqopds = 0.f;

    // Free track at the destination surface
    free_track_parameters<algebra_type> free_trk_1;
    free_trk_1.set_pos(r);
    free_trk_1.set_dir(t);
    free_trk_1.set_qop(free_trk_0.qop());

    // Path correction
    const free_matrix_t path_correction =
        jacobian_engine::template path_correction<destination_frame>(
            r, t, dtds, dqopds, trf_1);

    // Correction term for the path variation
    const free_matrix_t correction_term =
        matrix::identity<free_matrix_t>() + path_correction;

    // Free-to-bound jacobian at the destination surface
    const free_to_bound_matrix_t free_to_bound_jacobi =
        jacobian_engine::template free_to_bound_jacobian<destination_frame>(
            trf_1, free_trk_1);

    // Bound vector at the destination surface
    const bound_param_vector_t bound_vec_1 =
        detail::free_to_bound_vector<destination_frame>(trf_1, free_trk_1);

    // Full jacobian
    const bound_matrix_t full_jacobi = free_to_bound_jacobi * correction_term *
                                       transport_jacobi * bound_to_free_jacobi;

    // Update the covariance at the destination surface
    const bound_matrix_t bound_cov_1 =
        full_jacobi * bound_cov_0 * matrix::transpose(full_jacobi);

    bound_track_parameters<algebra_type> ret;
    ret.set_parameter_vector(bound_vec_1);
    ret.set_covariance(bound_cov_1);

    return ret;
  }

  intersection_t get_intersection(
      const std::array<intersection_t, 2>& inters) const {
    return inters[0];
  }

  intersection_t get_intersection(const intersection_t& inters) const {
    return inters;
  }
};

using TestTypes =
    ::testing::Types<annulus_type, rectangle_type, trapezoid_type, ring_type,
                     cylinder_type, straw_tube_type, drift_cell_type>;
TYPED_TEST_SUITE(detray_propagation_HelixCovarianceTransportValidation,
                 TestTypes, );

TYPED_TEST(detray_propagation_HelixCovarianceTransportValidation,
           one_loop_test) {
  // @NOTE: The test with high energy (>1 GeV) might fail due
  // to the numerical instability
  free_track_parameters<test_algebra> free_trk(
      {static_cast<scalar>(0.f), static_cast<scalar>(0.f),
       static_cast<scalar>(0.f)},
      static_cast<scalar>(0.f),
      {0.1f * unit<scalar>::GeV, static_cast<scalar>(0.f),
       static_cast<scalar>(0.f)},
      static_cast<scalar>(-1.f));

  detail::helix<test_algebra> reference_helix(free_trk, this->B);

  const std::size_t n_planes = 10u;
  std::vector<transform3> trfs =
      this->create_transforms(reference_helix, n_planes);
  ASSERT_EQ(trfs.size(), 10u);

  // Set the initial bound vector
  bound_parameters_vector<test_algebra> bound_vec_0 =
      detail::free_to_bound_vector<
          typename TestFixture::first_local_frame_type>(trfs[0], free_trk);

  // Set the initial bound covariance
  auto bound_cov_0 = matrix::zero<
      typename bound_track_parameters<test_algebra>::covariance_type>();
  getter::element(bound_cov_0, e_bound_loc0, e_bound_loc0) = 1.f;
  getter::element(bound_cov_0, e_bound_loc1, e_bound_loc1) = 1.f;
  getter::element(bound_cov_0, e_bound_phi, e_bound_phi) = 1.f;
  // Set theta error to zero, to suppress the loc1 (z) divergence
  getter::element(bound_cov_0, e_bound_theta, e_bound_theta) = 0.f;
  getter::element(bound_cov_0, e_bound_qoverp, e_bound_qoverp) = 1.f;
  getter::element(bound_cov_0, e_bound_time, e_bound_time) = 0.f;

  // Set bound track parameters
  bound_track_parameters<test_algebra> bound_params;
  bound_params.set_parameter_vector(bound_vec_0);
  bound_params.set_covariance(bound_cov_0);

  // Create masks
  const auto first_mask = std::get<rectangle_type>(this->masks);
  const auto test_mask = std::get<typename TestFixture::mask_type>(this->masks);

  // Total path length, just for testing purpose
  scalar total_path_length = 0.f;

  // Intersections for testing purporse
  std::vector<intersection_t> sfis;

  // Iterate over the planes until we reach the first plane (one loop)
  for (std::size_t i_p = 0u; i_p < n_planes; i_p++) {
    if (i_p == 0) {
      bound_params =
          this->propagate(bound_params, trfs[i_p], trfs[i_p + 1], first_mask,
                          test_mask, total_path_length, sfis);
    } else if (i_p > 0 && i_p < n_planes - 1) {
      bound_params =
          this->propagate(bound_params, trfs[i_p], trfs[i_p + 1], test_mask,
                          test_mask, total_path_length, sfis);

    } else if (i_p == n_planes - 1) {
      bound_params =
          this->propagate(bound_params, trfs[i_p], trfs[0], test_mask,
                          first_mask, total_path_length, sfis);
    }
  }

  // Check if the total path length is the expected value
  ASSERT_TRUE(total_path_length > 1e-3);
  ASSERT_NEAR(total_path_length,
              2.f * constant<scalar>::pi /
                  std::abs(reference_helix.qop() * reference_helix.B()),
              this->tolerance);

  ASSERT_EQ(sfis.size(), n_planes);
  for (std::size_t i = 0u; i < n_planes; i++) {
    EXPECT_TRUE(sfis[i].is_inside());
    EXPECT_TRUE(sfis[i].is_along());
  }

  // Check if the same vector is obtained after one loop
  for (unsigned int i = 0u; i < e_bound_size; i++) {
    EXPECT_NEAR(bound_vec_0[i], bound_params[i], this->tolerance)
        << "i: " << i << "\n"
        << std::setprecision(8) << bound_params;
  }

  // Check if the same covariance is obtained after one loop
  for (unsigned int i = 0u; i < e_bound_size; i++) {
    for (unsigned int j = 0u; j < e_bound_size; j++) {
      EXPECT_NEAR(getter::element(bound_cov_0, i, j),
                  getter::element(bound_params.covariance(), i, j),
                  this->tolerance)
          << "i: " << i << "\nj: " << j << "\n"
          << std::setprecision(8) << bound_params;
    }
  }
}

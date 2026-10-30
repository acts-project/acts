// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/core/concepts.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/tracking_surface.hpp"
#include "detray/navigation/caching_navigator.hpp"
#include "detray/navigation/intersection/helix_intersector.hpp"
#include "detray/propagator/actors/parameter_updater.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"

// Detray test include(s)
#include "detray/test/common/bfield.hpp"
#include "detray/test/common/build_telescope_detector.hpp"
#include "detray/test/common/track_generators.hpp"
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s).
#include <algorithm>
#include <iostream>
#include <optional>
#include <random>
#include <string>
#include <tuple>

using namespace detray;

// Type declarations
using test_algebra = test::algebra;
using scalar = dscalar<test_algebra>;
using vector3 = dvector3D<test_algebra>;
using transform3_type = dtransform3D<test_algebra>;
using bound_param_vector_type =
    bound_track_parameters<test_algebra>::parameter_vector_type;
using bound_covariance_type =
    bound_track_parameters<test_algebra>::covariance_type;
template <std::size_t ROWS, std::size_t COLS>
using matrix_type = dmatrix<test_algebra, ROWS, COLS>;

constexpr bool is_float = std::is_same_v<scalar, float>;

/// Test configuration
enum class field_type { e_const, e_inhom, e_inhom_material };
enum class surface_type { e_rect, e_wire };
using test_param = std::tuple<field_type, surface_type, navigation::direction>;

namespace {

// B-field Setup
// 1.996 T from averaged Bz strength of ODD field in the region of sqrt(x^2 +
// y^2) < 500 mm and abs(z) < 500 mm
const vector3 B_z{0.f, 0.f, 1.996f * unit<scalar>::T};

// Initial delta for numerical differentiaion
const std::array<scalar, 5u> h_sizes_rect{1e0f, 1e0f, 2e-2f, 1e-3f, 1e-3f};
const std::array<scalar, 5u> h_sizes_wire{1e0f, 1e0f, 2e-2f, 1e-3f, 1e-3f};

// Ridders' algorithm setup
constexpr const unsigned int Nt = 100u;
const std::array<scalar, 5u> safe{5.0f, 5.0f, 5.0f, 5.0f, 5.0f};
const std::array<scalar, 5u> con{1.1f, 1.1f, 1.1f, 1.1f, 1.1f};
constexpr const scalar big = std::numeric_limits<scalar>::max();

// Tolerances for the propagation
constexpr const scalar rk_tol_reference = 1e-4f * unit<scalar>::mm;
constexpr const scalar rk_tol_displaced = 1e-6f * unit<scalar>::mm;
constexpr const scalar on_surface_tol = 1e-3f * unit<scalar>::mm;
constexpr const scalar helix_tol = 1e-3f * unit<scalar>::mm;

// Number of tracks per configuration
constexpr const std::size_t n_tracks = 50u;

// Jacobian comparison: an element is accepted if it agrees relatively to
// jac_rtol, or if its absolute error moves the destination parameter i by less
// than resolution[i] under a one-sigma shift spread[j] of the departure
// parameter j. The absolute floor keeps elements that are zero or tiny in one
// of the two evaluations from being skipped or compared with an undefined
// relative error.
constexpr const scalar jac_rtol = is_float ? 1e-2f : 1e-3f;
const std::array<scalar, 5u> resolution{10.f * unit<scalar>::um,
                                        10.f * unit<scalar>::um, 1e-5f, 1e-5f,
                                        1e-5f /* relative to |qop| */};
const std::array<scalar, 5u> spread{
    1.f * unit<scalar>::mm, 1.f * unit<scalar>::mm, 1.f * unit<scalar>::mrad,
    1.f * unit<scalar>::mrad, 1e-2f /* relative to |qop| */};

// Covariance comparison, relative to sqrt(C_ii C_jj)
constexpr const scalar cov_rtol = is_float ? 2e-2f : 2e-3f;

std::random_device rd;
// For detector generation
std::mt19937_64 mt1(rd());
// For smearing initial parameter
std::mt19937_64 mt2(rd());

// Momentum range
constexpr const scalar min_mom = 0.5f * unit<scalar>::GeV;
constexpr const scalar max_mom = 100.f * unit<scalar>::GeV;

// Detector length generator
constexpr const scalar min_detector_length = 50.f * unit<scalar>::mm;
constexpr const scalar max_detector_length = 500.f * unit<scalar>::mm;
std::uniform_real_distribution rand_length(min_detector_length,
                                           max_detector_length);
constexpr const scalar envelope_size = 2000.f * unit<scalar>::mm;

// Mask size scaler
constexpr const scalar mask_scaler = 1.5f;

// Euler angles for the surface rotation
std::uniform_real_distribution<scalar> rand_alpha(0.f,
                                                  2.f * constant<scalar>::pi);
std::uniform_real_distribution<scalar> rand_cosbeta(constant<scalar>::inv_sqrt2,
                                                    1.f);
std::uniform_int_distribution rand_bool(0, 1);
std::uniform_real_distribution<scalar> rand_gamma(0.f,
                                                  2.f * constant<scalar>::pi);

// Particle setup = muon
pdg_particle ptc = muon<scalar>();

// Correlation factor in the range of [-10%, 10%]
constexpr scalar min_corr = -0.1f;
constexpr scalar max_corr = 0.1f;

// Values for sampling standard deviation
const std::array<scalar, 6u> stddevs_sampling{50.f * unit<scalar>::um,
                                              50.f * unit<scalar>::um,
                                              1.f * unit<scalar>::mrad,
                                              1.f * unit<scalar>::mrad,
                                              0.01f,
                                              1.f * unit<scalar>::ns};

// surface types
using rect_type = rectangle2D;
using wire_type = line_square;

}  // namespace

// Preprocess delta
void preprocess_delta(const unsigned int i, scalar& delta,
                      const bound_track_parameters<test_algebra> ref_param) {
  if (i == e_bound_theta) {
    const scalar rtheta = ref_param.theta();
    if (rtheta < constant<scalar>::pi_2) {
      delta = math::min(delta, rtheta);
    } else if (rtheta >= constant<scalar>::pi_2) {
      delta = math::min(delta, constant<scalar>::pi - rtheta);
    }
  }
}

// Containers for Ridders algorithm
// Page 231 of [Numerical Recipes] 3rd edition
struct ridders_derivative {
  std::array<std::array<std::array<scalar, Nt>, Nt>, 5u> Arr{};
  std::array<scalar, 5u> fac{};
  std::array<scalar, 5u> errt{};
  std::array<scalar, 5u> err{big, big, big, big, big};
  std::array<bool, 5u> complete{false, false, false, false, false};

  void initialize(const bound_param_vector_type& nvec1,
                  const bound_param_vector_type& nvec2, const scalar delta) {
    for (unsigned int j = 0; j < 5u; j++) {
      Arr[j][0][0] = (nvec1[j] - nvec2[j]) / (2.f * delta);
    }
  }

  void run(const bound_param_vector_type& nvec1,
           const bound_param_vector_type& nvec2, const scalar delta,
           const unsigned int p, const unsigned int i,
           bound_covariance_type& differentiated_jacobian) {
    for (unsigned int j = 0; j < 5u; j++) {
      Arr[j][0][p] = (nvec1[j] - nvec2[j]) / (2.f * delta);
    }

    const scalar con2 = con[i] * con[i];
    fac[i] = con2;

    for (unsigned int q = 1; q <= p; q++) {
      for (unsigned int j = 0; j < 5u; j++) {
        Arr[j][q][p] = (Arr[j][q - 1][p] * fac[i] - Arr[j][q - 1][p - 1]) /
                       (fac[i] - 1.0f);

        errt[j] = math::max(math::fabs(Arr[j][q][p] - Arr[j][q - 1][p]),
                            math::fabs(Arr[j][q][p] - Arr[j][q - 1][p - 1]));

        if (errt[j] <= err[j]) {
          if (complete[j] == false) {
            err[j] = errt[j];
            getter::element(differentiated_jacobian, j, i) = Arr[j][q][p];
            /*
            // Please leave this for debug
            if (j == e_bound_theta && i == e_bound_loc0) {
                std::clog << q << " " << p << " "
                          << getter::element(
                                 differentiated_jacobian, j, i)
                          << "  " << math::fabs(Arr[j][q][p])
                          << std::endl;
            }
            */
          }
        }
      }
      fac[i] = con2 * fac[i];
    }

    for (unsigned int j = 0; j < 5u; j++) {
      /*
      // Please leave this for debug
      if (j == e_bound_loc0 && i == e_bound_theta) {
          std::clog << getter::element(differentiated_jacobian, j, i)
                    << "  " << Arr[j][p][p] << "  "
                    << Arr[j][p - 1][p - 1] << "  "
                    << math::fabs(Arr[j][p][p] - Arr[j][p - 1][p - 1])
                    << "  " << safe[i] * err[j] << std::endl;
      }
      */
      if (math::fabs(Arr[j][p][p] - Arr[j][p - 1][p - 1]) >= safe[i] * err[j]) {
        complete[j] = true;
      }
    }
  }

  bool finished() { return (std::ranges::count(complete, false) == 0u); }
};

void wrap_angles(const bound_param_vector_type& ref_vector,
                 bound_param_vector_type& target_vector) {
  const scalar rphi = ref_vector.phi();
  const scalar tphi = target_vector.phi();
  scalar new_tphi = tphi;

  if (rphi >= constant<scalar>::pi_2) {
    if (tphi < 0.f) {
      new_tphi += 2.f * constant<scalar>::pi;
    }
  } else if (rphi < -constant<scalar>::pi_2) {
    if (tphi >= 0.f) {
      new_tphi -= 2.f * constant<scalar>::pi;
    }
  }

  target_vector.set_phi(new_tphi);
}

// Get random initial covariance
bound_covariance_type get_random_initial_covariance(const scalar ini_qop) {
  // Initial covariance matrix for smearing
  auto ini_cov = matrix::zero<bound_covariance_type>();

  // Random correction factor
  std::uniform_real_distribution<scalar> rand_corr(min_corr, max_corr);

  // Distribution for sampling standard deviations
  std::normal_distribution<scalar> rand_l0(0.f * unit<scalar>::um,
                                           stddevs_sampling[0u]);
  std::normal_distribution<scalar> rand_l1(0.f * unit<scalar>::um,
                                           stddevs_sampling[1u]);
  std::normal_distribution<scalar> rand_phi(0.f * unit<scalar>::mrad,
                                            stddevs_sampling[2u]);
  std::normal_distribution<scalar> rand_theta(0.f * unit<scalar>::mrad,
                                              stddevs_sampling[3u]);
  std::normal_distribution<scalar> rand_qop(0.f,
                                            stddevs_sampling[4u] * ini_qop);
  std::normal_distribution<scalar> rand_time(0.f * unit<scalar>::ns,
                                             stddevs_sampling[5u]);

  std::array<scalar, 6u> stddevs{};
  stddevs[0] = rand_l0(mt2);
  stddevs[1] = rand_l1(mt2);
  stddevs[2] = rand_phi(mt2);
  stddevs[3] = rand_theta(mt2);
  stddevs[4] = rand_qop(mt2);
  stddevs[5] = rand_time(mt2);

  for (unsigned int i = 0u; i < 6u; i++) {
    for (unsigned int j = 0u; j < 6u; j++) {
      if (i == j) {
        getter::element(ini_cov, i, i) = stddevs[i] * stddevs[i];
      } else if (i > j) {
        getter::element(ini_cov, i, j) =
            std::abs(stddevs[i]) * std::abs(stddevs[j]) * rand_corr(mt2);
        getter::element(ini_cov, j, i) = getter::element(ini_cov, i, j);
      }
    }
  }

  return ini_cov;
}

template <concepts::detector detector_t,
          typename detector_t::metadata::mask_id mask_id>
std::pair<euler_rotation<test_algebra>, std::array<scalar, 3u>> tilt_surface(
    detector_t& det, const unsigned int sf_id, const vector3& helix_dir,
    const scalar alpha, const scalar beta, const scalar gamma) {
  const auto& sf = det.surface(sf_id);
  const auto& trf_link = sf.transform();
  auto& trf = det.transform_store().at(trf_link);

  euler_rotation<test_algebra> euler;
  euler.alpha = alpha;

  if (sf_id == 1u) {
    euler.beta = beta;
    euler.gamma = gamma;
  }

  // Helix direction
  euler.z = helix_dir;

  if constexpr (mask_id == detector_t::masks::id::e_rectangle2D) {
    // ubasis = trf.x() for bound frame
    euler.x = trf.x();
  } else if (mask_id == detector_t::masks::id::e_drift_cell) {
    // ubasis = vxt/|vxt| where v is trf.z() and t is helix_dir
    EXPECT_NEAR(vector::dot(trf.z(), helix_dir), 0.f, 1e-6f);
    euler.x = vector::normalize(vector::cross(trf.z(), helix_dir));
  }

  EXPECT_NEAR(vector::dot(euler.x, euler.z), 0.f, 1e-6f);

  auto [local_x, local_z] = euler();

  if constexpr (mask_id == detector_t::masks::id::e_drift_cell) {
    // local_z should be the line direction
    local_z = vector::cross(local_z, local_x);
  }

  // At the moment, we do not shift the surfaces
  scalar x_shift{0.f};
  scalar y_shift{0.f};
  scalar z_shift{0.f};

  // Translation vector
  vector3 translation = trf.translation();

  vector3 displacement({x_shift, y_shift, z_shift});
  translation = trf.translation() + displacement;
  transform3_type new_trf(translation, local_z, local_x, true);
  // @todo : Remove and use different context
  detray::detail::set_transform(det, new_trf, trf_link);

  return {euler, {x_shift, y_shift, z_shift}};
}

template <concepts::algebra algebra_t>
struct bound_getter : public base_actor {
  // Track types
  using bound_track_parameters_type = bound_track_parameters<algebra_t>;
  using free_track_parameters_type = free_track_parameters<algebra_t>;

  struct state {
    unsigned int m_departure_idx{0u};
    unsigned int m_destination_idx{1u};
    bool m_reached_departure{false};
    bool m_reached_destination{false};
    scalar m_min_path_length{detray::detail::invalid_value<scalar>()};
    scalar m_path_length{0.f};
    scalar m_abs_path_length{0.f};
    bound_track_parameters_type m_param_departure{};
    bound_track_parameters_type m_param_destination{};
    scalar m_avg_step_size{0.f};
    std::size_t step_count{0u};
    std::size_t track_ID{0u};
  };

  template <typename propagator_state_t>
  DETRAY_HOST_DEVICE void operator()(
      state& actor_state, propagator_state_t& propagation,
      const detray::actor::parameter_transporter_result<algebra_t>& res) const {
    auto& navigation = propagation.navigation();
    auto& stepping = propagation.stepping();

    actor_state.step_count++;

    const scalar N = static_cast<scalar>(actor_state.step_count);

    actor_state.m_avg_step_size =
        ((N - 1.f) * actor_state.m_avg_step_size + stepping.step_size()) / N;

    // Warning for too many step counts
    if (actor_state.step_count > 1000000) {
      std::clog << "Too many step counts!" << std::endl;
      std::clog << "Track ID: " << actor_state.track_ID << std::endl;
      std::clog << "Path length: " << actor_state.m_path_length << std::endl;
      std::clog << "PhiI: " << actor_state.m_param_departure.phi() << std::endl;
      std::clog << "ThetaI: " << actor_state.m_param_departure.theta()
                << std::endl;
      std::clog << "QopI: " << actor_state.m_param_departure.qop() << std::endl;
      navigation.exit();
      propagation.heartbeat(false);
    }

    if ((navigation.is_on_sensitive() || navigation.is_on_passive()) &&
        navigation.geometry_identifier().index() ==
            actor_state.m_departure_idx &&
        !actor_state.m_reached_departure) {
      actor_state.m_param_departure = res.destination_params();
      actor_state.m_reached_departure = true;
    }
    // Get the bound track parameters at the destination surface. The full
    // Jacobian is not recomputed here: by the time an observer runs, the
    // parameter updater has already reset the transport Jacobian. It is
    // instead accumulated by the updater itself through its Jacobian pointer.
    else if ((navigation.is_on_sensitive() || navigation.is_on_passive()) &&
             navigation.geometry_identifier().index() ==
                 actor_state.m_destination_idx &&
             actor_state.m_reached_departure) {
      actor_state.m_path_length = stepping.path_length();
      actor_state.m_abs_path_length = stepping.abs_path_length();
      actor_state.m_param_destination = res.destination_params();
      actor_state.m_reached_destination = true;

      // Stop navigation if the destination surface is found
      navigation.exit();
      propagation.heartbeat(false);
    }

    if (math::fabs(stepping.path_length()) > actor_state.m_min_path_length) {
      propagation.navigation().set_no_trust();
    }

    return;
  }
};

/// Numerically integrate the jacobian
template <typename propagator_t, typename field_t>
bound_getter<test_algebra>::state evaluate_bound_param(
    const std::size_t trk_count, const scalar detector_length,
    const bound_track_parameters<test_algebra>& initial_param,
    const typename propagator_t::detector_type& det, const field_t& field,
    const scalar overstep_tolerance, const scalar path_tolerance,
    const scalar rk_tolerance, const navigation::direction dir,
    bool use_field_gradient, bool do_covariance_transport, bool do_inspect,
    bound_covariance_type* full_jacobian = nullptr) {
  // Propagator is built from the stepper and navigator
  propagation::config cfg{};
  cfg.navigation.intersection.overstep_tolerance =
      static_cast<float>(overstep_tolerance);
  cfg.navigation.intersection.path_tolerance =
      static_cast<float>(path_tolerance);
  cfg.navigation.estimate_scattering_noise = false;
  cfg.stepping.rk_error_tol = static_cast<float>(rk_tolerance);
  cfg.stepping.use_eloss_gradient = true;
  cfg.stepping.use_field_gradient = use_field_gradient;
  cfg.stepping.do_covariance_transport = do_covariance_transport;
  propagator_t p(cfg);

  // Actor states. The full Jacobian is accumulated by the updater through the
  // same pointer mechanism that downstream clients use.
  actor::parameter_updater_state<test_algebra> updater_state{cfg, initial_param,
                                                             full_jacobian};
  bound_getter<test_algebra>::state bound_getter_state{};
  bound_getter_state.track_ID = trk_count;
  bound_getter_state.m_min_path_length = detector_length * 0.75f;
  bound_getter_state.m_departure_idx =
      (dir == navigation::direction::e_forward) ? 0u : 1u;
  bound_getter_state.m_destination_idx =
      1u - bound_getter_state.m_departure_idx;
  auto actor_states = detray::tie(updater_state, bound_getter_state);

  // Init propagator states for the reference track
  typename propagator_t::stepper_type::magnetic_field_type bfield_view(field);
  typename propagator_t::state state(initial_param, bfield_view, det);

  // Run the propagation for the reference track
  state.set_particle(ptc);
  state.debug(do_inspect);
  state.navigation().set_direction(dir);

  p.propagate(state, actor_states);

  return bound_getter_state;
}

template <typename propagator_t, typename field_t>
bound_param_vector_type get_displaced_bound_vector(
    const std::size_t trk_count,
    const bound_track_parameters<test_algebra>& ref_param,
    const typename propagator_t::detector_type& det,
    const scalar detector_length, const field_t& field,
    const scalar overstep_tolerance, const scalar path_tolerance,
    const scalar rk_tolerance, const navigation::direction dir,
    const unsigned int target_index, const scalar displacement, bool& ok) {
  propagation::config cfg{};
  cfg.navigation.intersection.overstep_tolerance =
      static_cast<float>(overstep_tolerance);
  cfg.navigation.intersection.path_tolerance =
      static_cast<float>(path_tolerance);
  cfg.navigation.estimate_scattering_noise = false;
  cfg.stepping.rk_error_tol = static_cast<float>(rk_tolerance);
  cfg.stepping.use_eloss_gradient = true;
  cfg.stepping.do_covariance_transport = false;

  // Propagator is built from the stepper and navigator
  propagator_t p(cfg);

  bound_track_parameters<test_algebra> dparam = ref_param;
  dparam[target_index] += displacement;

  typename propagator_t::stepper_type::magnetic_field_type bfield_view(field);
  typename propagator_t::state dstate(dparam, bfield_view, det);

  // Actor states
  actor::parameter_updater_state<test_algebra> updater_state{cfg, dparam};
  bound_getter<test_algebra>::state bound_getter_state{};
  bound_getter_state.track_ID = trk_count;
  bound_getter_state.m_min_path_length = detector_length * 0.75f;
  bound_getter_state.m_departure_idx =
      (dir == navigation::direction::e_forward) ? 0u : 1u;
  bound_getter_state.m_destination_idx =
      1u - bound_getter_state.m_departure_idx;

  auto actor_states = detray::tie(updater_state, bound_getter_state);
  dstate.set_particle(ptc);
  dstate.navigation().set_direction(dir);

  p.propagate(dstate, actor_states);

  ok = ok && bound_getter_state.m_reached_destination;
  bound_param_vector_type new_vec = bound_getter_state.m_param_destination;

  // phi needs to be wrapped w.r.t. phi of the reference vector
  wrap_angles(ref_param, new_vec);

  return new_vec;
}

/// Numerically evaluate the jacobian
template <typename propagator_t, typename field_t>
bound_track_parameters<test_algebra>::covariance_type directly_differentiate(
    const std::size_t trk_count,
    const bound_track_parameters<test_algebra>& ref_param,
    const typename propagator_t::detector_type& det,
    const scalar detector_length, const field_t& field,
    const scalar overstep_tolerance, const scalar path_tolerance,
    const scalar rk_tolerance, const navigation::direction dir,
    const std::array<scalar, 5u> hs,
    std::array<unsigned int, 5u>& num_iterations,
    std::array<bool, 25>& convergence, bool& ok) {
  // Return Jacobian
  bound_covariance_type differentiated_jacobian =
      matrix::identity<bound_covariance_type>();

  for (unsigned int i = 0u; i < 5u; i++) {
    scalar delta = hs[i];
    preprocess_delta(i, delta, ref_param);

    const auto vec1 = get_displaced_bound_vector<propagator_t, field_t>(
        trk_count, ref_param, det, detector_length, field, overstep_tolerance,
        path_tolerance, rk_tolerance, dir, i, 1.f * delta, ok);
    const auto vec2 = get_displaced_bound_vector<propagator_t, field_t>(
        trk_count, ref_param, det, detector_length, field, overstep_tolerance,
        path_tolerance, rk_tolerance, dir, i, -1.f * delta, ok);

    ridders_derivative ridder;
    ridder.initialize(vec1, vec2, delta);

    for (unsigned int p = 1u; p < Nt; p++) {
      delta /= con[i];

      const auto nvec1 = get_displaced_bound_vector<propagator_t, field_t>(
          trk_count, ref_param, det, detector_length, field, overstep_tolerance,
          path_tolerance, rk_tolerance, dir, i, 1.f * delta, ok);
      const auto nvec2 = get_displaced_bound_vector<propagator_t, field_t>(
          trk_count, ref_param, det, detector_length, field, overstep_tolerance,
          path_tolerance, rk_tolerance, dir, i, -1.f * delta, ok);

      ridder.run(nvec1, nvec2, delta, p, i, differentiated_jacobian);

      if (ridder.finished()) {
        num_iterations[i] = p;
        break;
      }
    }

    // Row-major
    for (std::size_t j = 0u; j < 5u; j++) {
      convergence[i + j * 5] = ridder.complete[j];
    }
  }

  return differentiated_jacobian;
}

template <concepts::detector detector_t,
          typename detector_t::metadata::mask_id mid>
bound_track_parameters<test_algebra> get_initial_parameter(
    const detector_t& det, const free_track_parameters<test_algebra>& vertex,
    const vector3& field, const scalar helix_tolerance) {
  // Helix from the vertex
  detail::helix<test_algebra> hlx(vertex, field);

  const auto& departure_sf = det.surface(0u);
  const auto& trf_link = departure_sf.transform();
  const auto& departure_trf = det.transform_store().at(trf_link);
  const auto& mask_link = departure_sf.mask();
  const auto& departure_mask =
      det.mask_store().template get<mid>().at(mask_link.index());

  using mask_t = types::get<typename detector_t::masks, mid>;
  helix_intersector<typename mask_t::shape, test_algebra> hlx_is{};
  hlx_is.run_rtsafe = false;
  hlx_is.convergence_tolerance = helix_tolerance;
  auto sfi = hlx_is(hlx, departure_sf, departure_mask, departure_trf, 0.f);
  EXPECT_TRUE(sfi.is_inside())
      << " Initial surface not found" << std::endl
      << " log10(Helix tolerance): " << math::log10(helix_tolerance)
      << " Phi: " << vector::phi(vertex.dir())
      << " Theta: " << vector::theta(vertex.dir())
      << " Mom [GeV/c]: " << vertex.p(ptc.charge()) << std::endl
      << sfi;

  const auto path_length = sfi.path();
  // As we don't rotate or shift the initial surface anymore, the path_length
  // should be 0
  EXPECT_FLOAT_EQ(static_cast<float>(path_length), 0.f);

  const auto pos = hlx(path_length);
  const auto dir = hlx.dir(path_length);

  const free_track_parameters<test_algebra> free_par(pos, 0, dir, hlx.qop());

  const auto bound_vec =
      tracking_surface{det, departure_sf}.free_to_bound_vector({}, free_par);

  bound_track_parameters<test_algebra> ret;
  ret.set_surface_link(departure_sf.identifier());
  ret.set_parameter_vector(bound_vec);

  return ret;
}

/// Tolerance for element (i, j) of the Jacobian
scalar jacobian_tolerance(const unsigned int i, const unsigned int j,
                          const scalar numerical, const scalar qop) {
  scalar res_i = resolution[i];
  scalar spread_j = spread[j];
  if (i == e_bound_qoverp) {
    res_i *= math::fabs(qop);
  }
  if (j == e_bound_qoverp) {
    spread_j *= math::fabs(qop);
  }
  return jac_rtol * math::fabs(numerical) + res_i / spread_j;
}

const char* param_name(const unsigned int i) {
  constexpr std::array<const char*, 6u> names{"l0",    "l1",  "phi",
                                              "theta", "qop", "t"};
  return names[i];
}

/// Validate the bound-to-bound Jacobian and the covariance transport of the
/// propagator against central finite differences (Ridders' method) of the
/// propagation map between the two surfaces of a randomly tilted telescope.
class JacobianValidation : public ::testing::TestWithParam<test_param> {};

TEST_P(JacobianValidation, propagator) {
  const auto [field_cfg, surface_cfg, dir] = GetParam();

  // The inhomogeneous field map is read from the environment, as in the
  // other detray tests
  if (field_cfg != field_type::e_const &&
      std::getenv("DETRAY_BFIELD_FILE") == nullptr) {
    GTEST_SKIP() << "DETRAY_BFIELD_FILE is not set";
  }

  vecmem::host_memory_resource host_mr;

  // Deterministic per configuration
  const std::size_t mc_seed =
      42u + 8u * static_cast<std::size_t>(field_cfg) +
      2u * static_cast<std::size_t>(surface_cfg) +
      static_cast<std::size_t>(dir == navigation::direction::e_backward);
  mt1.seed(mc_seed);
  mt2.seed(mc_seed);

  // Fields
  const auto const_field = create_const_field<scalar>(B_z);
  std::optional<bfield::inhom_field_t<scalar>> inhom_field;
  if (field_cfg != field_type::e_const) {
    inhom_field = create_inhom_field<scalar>();
  }

  // Volume material
  material<scalar> volume_mat = vacuum<scalar>{};
  if (field_cfg == field_type::e_inhom_material) {
    volume_mat = cesium_iodide_with_ded<scalar>{};
  }

  // Track generator: random direction and momentum from the origin
  using uniform_gen_t =
      detail::random_numbers<scalar, std::uniform_real_distribution<scalar>>;
  using trk_generator_t =
      random_track_generator<free_track_parameters<test_algebra>,
                             uniform_gen_t>;
  trk_generator_t::configuration trk_gen_cfg{};
  trk_gen_cfg.seed(mc_seed);
  trk_gen_cfg.n_tracks(n_tracks);
  trk_gen_cfg.phi_range(-constant<scalar>::pi, constant<scalar>::pi);
  trk_gen_cfg.theta_range(0.f, constant<scalar>::pi);
  trk_gen_cfg.mom_range(min_mom, max_mom);
  trk_gen_cfg.origin(0.f, 0.f, 0.f);
  trk_gen_cfg.origin_stddev(0.f, 0.f, 0.f);

  // Runs one configuration for a given surface shape
  auto run = [&]<typename shape_t, typename field_t>(const field_t& field) {
    using detector_t =
        host::detector<telescope_metadata<test_algebra, shape_t>>;
    using navigator_t = caching_navigator<detector_t>;
    using stepper_t =
        rk_stepper<typename field_t::view_t, test_algebra,
                   constrained_step<scalar>, stepper_default_policy<scalar>>;
    using actor_chain_t = actor_chain<
        actor::parameter_updater<test_algebra, bound_getter<test_algebra>>>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

    constexpr bool is_wire = std::is_same_v<shape_t, wire_type>;
    constexpr auto mask_id = is_wire ? detector_t::masks::id::e_drift_cell
                                     : detector_t::masks::id::e_rectangle2D;
    const std::array<scalar, 5u>& hs = is_wire ? h_sizes_wire : h_sizes_rect;

    std::size_t n_checked_elements{0u};
    std::size_t n_unconverged_elements{0u};
    std::size_t track_count{0u};

    for (const auto track : trk_generator_t{trk_gen_cfg}) {
      const std::size_t track_id = track_count++;
      SCOPED_TRACE("track " + std::to_string(track_id));

      // Pilot track
      detail::helix<test_algebra> helix_bz(track, B_z);

      // Make a telescope geometry around the pilot track
      const scalar detector_length = rand_length(mt1);

      mask<shape_t, test_algebra> m{0u, detector_length * mask_scaler,
                                    detector_length * mask_scaler};

      // Adjust overstep tolerance
      const scalar overstep_tol =
          -math::min(1000.f * unit<scalar>::mm, detector_length * 0.5f);

      tel_det_config<test_algebra, shape_t, detail::helix> tel_cfg{m, helix_bz};
      tel_cfg.envelope(envelope_size);
      tel_cfg.module_material(vacuum<scalar>{});
      tel_cfg.mat_thickness(0.f);
      tel_cfg.n_surfaces(2u);
      tel_cfg.length(detector_length);
      tel_cfg.volume_material(volume_mat);
      tel_cfg.do_check(false);

      auto alphaI = rand_alpha(mt1);
      auto alphaF = rand_alpha(mt1);
      auto betaF = math::acos(rand_cosbeta(mt1));
      if (rand_bool(mt1) == 0) {
        betaF = -betaF;
      }
      auto gammaF = rand_gamma(mt1);

      auto [det, names] =
          build_telescope_detector<test_algebra>(host_mr, tel_cfg);
      tilt_surface<detector_t, mask_id>(det, 0u, helix_bz.dir(0.f), alphaI, 0.f,
                                        0.f);
      tilt_surface<detector_t, mask_id>(det, 1u, helix_bz.dir(detector_length),
                                        alphaF, betaF, gammaF);

      // Start parameters on surface 0; for the backward direction, first go
      // forward to surface 1 and start from there
      bound_track_parameters<test_algebra> start =
          get_initial_parameter<detector_t, mask_id>(det, track, B_z,
                                                     helix_tol);

      if (dir == navigation::direction::e_backward) {
        const auto fw = evaluate_bound_param<propagator_t, field_t>(
            track_id, detector_length, start, det, field, overstep_tol,
            on_surface_tol, rk_tol_displaced, navigation::direction::e_forward,
            false, false, false);
        ASSERT_TRUE(fw.m_reached_destination)
            << "Forward pass to surface 1 failed";
        start = fw.m_param_destination;
      }
      start.set_covariance(get_random_initial_covariance(start.qop()));

      // Reference propagation with covariance transport and full Jacobian
      bound_covariance_type J = matrix::identity<bound_covariance_type>();
      const auto ref = evaluate_bound_param<propagator_t, field_t>(
          track_id, detector_length, start, det, field, overstep_tol,
          on_surface_tol, rk_tol_reference, dir, false, true, false, &J);
      ASSERT_TRUE(ref.m_reached_destination)
          << "Reference propagation did not reach the destination";
      ASSERT_GE(math::fabs(ref.m_path_length), 0.5f * detector_length);
      ASSERT_LE(math::fabs(ref.m_path_length), 1.5f * detector_length);

      // Finite differences of the same map, without covariance transport
      std::array<unsigned int, 5u> num_iterations{};
      std::array<bool, 25u> convergence{};
      bool fd_ok{true};
      const bound_covariance_type J_fd =
          directly_differentiate<propagator_t, field_t>(
              track_id, ref.m_param_departure, det, detector_length, field,
              overstep_tol, on_surface_tol, rk_tol_displaced, dir, hs,
              num_iterations, convergence, fd_ok);
      ASSERT_TRUE(fd_ok)
          << "A displaced propagation did not reach the destination";

      const scalar qop = ref.m_param_departure.qop();

      // Time is not propagated: its row and column must be trivial
      for (unsigned int k = 0u; k < 5u; k++) {
        EXPECT_EQ(getter::element(J, k, e_bound_time), 0.f) << param_name(k);
        EXPECT_EQ(getter::element(J, e_bound_time, k), 0.f) << param_name(k);
      }
      EXPECT_EQ(getter::element(J, e_bound_time, e_bound_time), 1.f);

      // All 25 elements of the 5x5 block (convergence is stored row-major as
      // convergence[i + j * 5] for d(param_j)/d(param_i))
      for (unsigned int j = 0u; j < 5u; j++) {
        for (unsigned int i = 0u; i < 5u; i++) {
          if (!convergence[i + j * 5u]) {
            n_unconverged_elements++;
            continue;
          }
          n_checked_elements++;
          const scalar num = getter::element(J_fd, j, i);
          const scalar ana = getter::element(J, j, i);
          EXPECT_NEAR(ana, num, jacobian_tolerance(j, i, num, qop))
              << "d" << param_name(j) << "/d" << param_name(i) << " (path "
              << ref.m_path_length << " mm, qop " << qop << ")";
        }
      }

      // Transported covariance against J_fd C J_fd^T. The finite-difference
      // Jacobian is 5x5; the time row and column are the identity.
      bound_covariance_type J6 = J_fd;
      for (unsigned int k = 0u; k < 5u; k++) {
        getter::element(J6, k, e_bound_time) = 0.f;
        getter::element(J6, e_bound_time, k) = 0.f;
      }
      getter::element(J6, e_bound_time, e_bound_time) = 1.f;

      const bound_covariance_type& C0 = start.covariance();
      const bound_covariance_type C_fd = J6 * C0 * matrix::transpose(J6);
      const bound_covariance_type& C = ref.m_param_destination.covariance();

      for (unsigned int i = 0u; i < e_bound_size; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
          const scalar scale = math::sqrt(getter::element(C_fd, i, i) *
                                          getter::element(C_fd, j, j));
          EXPECT_NEAR(getter::element(C, i, j), getter::element(C_fd, i, j),
                      cov_rtol * scale)
              << "C(" << param_name(i) << "," << param_name(j) << ")";
        }
      }

      // The transported covariance must be symmetric to rounding
      for (unsigned int i = 0u; i < e_bound_size; i++) {
        for (unsigned int j = 0u; j < i; j++) {
          const scalar scale =
              math::sqrt(getter::element(C, i, i) * getter::element(C, j, j));
          EXPECT_NEAR(getter::element(C, i, j), getter::element(C, j, i),
                      1e-5f * scale)
              << "asymmetry in C(" << param_name(i) << "," << param_name(j)
              << ")";
        }
      }
    }

    // The differentiation must converge for essentially all elements
    EXPECT_EQ(track_count, n_tracks);
    EXPECT_LE(n_unconverged_elements,
              (n_checked_elements + n_unconverged_elements) / 50u)
        << "Too many non-converged finite differences";
  };

  auto dispatch_field = [&]<typename shape_t>() {
    if (field_cfg == field_type::e_const) {
      run.template operator()<shape_t>(const_field);
    } else {
      run.template operator()<shape_t>(*inhom_field);
    }
  };

  if (surface_cfg == surface_type::e_rect) {
    dispatch_field.template operator()<rect_type>();
  } else {
    dispatch_field.template operator()<wire_type>();
  }
}

std::string test_name(const ::testing::TestParamInfo<test_param>& info) {
  const field_type f = std::get<0>(info.param);
  const surface_type sf = std::get<1>(info.param);
  const navigation::direction d = std::get<2>(info.param);
  std::string name = (f == field_type::e_const)   ? "const_field"
                     : (f == field_type::e_inhom) ? "inhom_field"
                                                  : "inhom_field_material";
  name += (sf == surface_type::e_rect) ? "_rect" : "_wire";
  name += (d == navigation::direction::e_forward) ? "_forward" : "_backward";
  return name;
}

INSTANTIATE_TEST_SUITE_P(
    detray_propagator, JacobianValidation,
    ::testing::Combine(
        ::testing::Values(field_type::e_const, field_type::e_inhom,
                          field_type::e_inhom_material),
        ::testing::Values(surface_type::e_rect, surface_type::e_wire),
        ::testing::Values(navigation::direction::e_forward,
                          navigation::direction::e_backward)),
    test_name);

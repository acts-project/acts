// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
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
#include "detray/test/utils/inspectors.hpp"
#include "detray/test/utils/statistics.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// Boost
#include "detray/options/boost_program_options.hpp"

// System include(s).
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace po = boost::program_options;
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
  std::array<std::array<std::array<scalar, Nt>, Nt>, 5u> Arr;
  std::array<scalar, 5u> fac;
  std::array<scalar, 5u> errt;
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
        fac[i] = con2 * fac[i];

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

scalar get_relative_difference(scalar ref_val, scalar num_val) {
  scalar rel_diff{0.f};

  // If the evaluated jacovian or numerical differentiation is too small set the
  // relative difference to zero
  if (ref_val == 0.f || num_val == 0.f) {
    rel_diff = 0.f;
  } else {
    rel_diff = std::abs(ref_val - num_val) / std::abs(num_val);
  }

  return rel_diff;
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

  std::array<scalar, 6u> stddevs;
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

// Input covariance should be the diagonal matrix
auto get_smeared_bound_vector(const bound_covariance_type& ini_cov,
                              const bound_param_vector_type& ini_vec) {
  using bound_vector_type = bound_vector<test_algebra>;

  // Do the Cholesky Decomposition
  const bound_covariance_type L = matrix::cholesky_decomposition(ini_cov);

  // Vecor with random elements from a normal distribution
  auto k = matrix::zero<bound_vector_type>();
  std::normal_distribution<scalar> normal_dist(0.f, 1.f);
  for (unsigned int i = 0u; i < 5u; i++) {
    // Smear the value
    getter::element(k, i, 0) = normal_dist(mt2);
  }

  const bound_vector_type new_vec = ini_vec.vector() + L * k;

  return bound_param_vector_type{new_vec};
}

template <typename detector_t, typename detector_t::metadata::mask_id mask_id>
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
    scalar m_min_path_length;
    scalar m_path_length;
    scalar m_abs_path_length;
    bound_track_parameters_type m_param_departure;
    bound_track_parameters_type m_param_destination;
    typename bound_track_parameters_type::covariance_type m_jacobi;
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
        navigation.geometry_identifier().index() == 0u) {
      actor_state.m_param_departure = res.destination_params();
    }
    // Get the bound track parameters and jacobian at the destination
    // surface
    else if ((navigation.is_on_sensitive() || navigation.is_on_passive()) &&
             navigation.geometry_identifier().index() == 1u) {
      actor_state.m_path_length = stepping.path_length();
      actor_state.m_abs_path_length = stepping.abs_path_length();
      actor_state.m_param_destination = res.destination_params();
      actor_state.m_jacobi =
          actor::parameter_transporter<algebra_t>().get_full_jacobian(
              propagation, actor_state.m_param_departure);

      // Stop navigation if the destination surface is found
      navigation.exit();
      propagation.heartbeat(false);
    }

    if (stepping.path_length() > actor_state.m_min_path_length) {
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
    const scalar rk_tolerance, const scalar constraint_step,
    bool use_field_gradient, bool do_covariance_transport, bool do_inspect) {
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

  // Actor states
  actor::parameter_updater_state<test_algebra> updater_state{cfg,
                                                             initial_param};
  bound_getter<test_algebra>::state bound_getter_state{};
  bound_getter_state.track_ID = trk_count;
  bound_getter_state.m_min_path_length = detector_length * 0.75f;
  auto actor_states = detray::tie(updater_state, bound_getter_state);

  // Init propagator states for the reference track
  typename propagator_t::stepper_type::magnetic_field_type bfield_view(field);
  typename propagator_t::state state(initial_param, bfield_view, det);

  // Run the propagation for the reference track
  state.set_particle(ptc);
  state.debug(do_inspect);
  state.stepping()
      .template set_constraint<detray::step::constraint::e_accuracy>(
          static_cast<float>(constraint_step));

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
    const scalar rk_tolerance, const scalar constraint_step,
    const unsigned int target_index, const scalar displacement) {
  propagation::config cfg{};
  cfg.navigation.intersection.overstep_tolerance =
      static_cast<float>(overstep_tolerance);
  cfg.navigation.intersection.path_tolerance =
      static_cast<float>(path_tolerance);
  cfg.navigation.estimate_scattering_noise = false;
  cfg.stepping.rk_error_tol = static_cast<float>(rk_tolerance);
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

  auto actor_states = detray::tie(updater_state, bound_getter_state);
  dstate.set_particle(ptc);
  dstate.stepping()
      .template set_constraint<detray::step::constraint::e_accuracy>(
          constraint_step);

  p.propagate(dstate, actor_states);

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
    const scalar rk_tolerance, const scalar constraint_step,
    const std::array<scalar, 5u> hs,
    std::array<unsigned int, 5u>& num_iterations,
    std::array<bool, 25>& convergence) {
  // Return Jacobian
  bound_covariance_type differentiated_jacobian;

  for (unsigned int i = 0u; i < 5u; i++) {
    scalar delta = hs[i];
    preprocess_delta(i, delta, ref_param);

    const auto vec1 = get_displaced_bound_vector<propagator_t, field_t>(
        trk_count, ref_param, det, detector_length, field, overstep_tolerance,
        path_tolerance, rk_tolerance, constraint_step, i, 1.f * delta);
    const auto vec2 = get_displaced_bound_vector<propagator_t, field_t>(
        trk_count, ref_param, det, detector_length, field, overstep_tolerance,
        path_tolerance, rk_tolerance, constraint_step, i, -1.f * delta);

    ridders_derivative ridder;
    ridder.initialize(vec1, vec2, delta);

    for (unsigned int p = 1u; p < Nt; p++) {
      delta /= con[i];

      const auto nvec1 = get_displaced_bound_vector<propagator_t, field_t>(
          trk_count, ref_param, det, detector_length, field, overstep_tolerance,
          path_tolerance, rk_tolerance, constraint_step, i, 1.f * delta);
      const auto nvec2 = get_displaced_bound_vector<propagator_t, field_t>(
          trk_count, ref_param, det, detector_length, field, overstep_tolerance,
          path_tolerance, rk_tolerance, constraint_step, i, -1.f * delta);

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

template <typename detector_t, typename detector_t::metadata::mask_id mid>
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
  ret.set_surface_link(geometry::identifier{0u});
  ret.set_parameter_vector(bound_vec);

  return ret;
}

template <typename propagator_t, typename field_t>
void evaluate_jacobian_difference(
    const std::size_t trk_count, const std::array<scalar, 3u>& euler_angles_I,
    const std::array<scalar, 3u>& euler_angles_F,
    const typename propagator_t::detector_type& det,
    const scalar detector_length,
    const bound_track_parameters<test_algebra>& track, const field_t& field,
    const scalar overstep_tolerance, const scalar path_tolerance,
    const scalar rk_tolerance, const scalar rk_tolerance_dis,
    const scalar constraint_step, const std::array<scalar, 5u>& hs,
    std::ofstream& file, scalar& ref_rel_diff, bool use_field_gradient,
    bool do_inspect, const bool use_precal_values = false,
    [[maybe_unused]] bound_covariance_type precal_diff_jacobi = {},
    [[maybe_unused]] std::array<unsigned int, 5u> precal_num_iterations = {},
    [[maybe_unused]] std::array<bool, 25u> precal_convergence = {}) {
  const auto phi0 = track.phi();
  const auto theta0 = track.theta();
  (void)phi0;
  (void)theta0;

  auto bound_getter = evaluate_bound_param<propagator_t, field_t>(
      trk_count, detector_length, track, det, field, overstep_tolerance,
      path_tolerance, rk_tolerance, constraint_step, use_field_gradient, true,
      do_inspect);

  const auto reference_param = bound_getter.m_param_departure;
  const auto final_param = bound_getter.m_param_destination;

  // Sanity check
  ASSERT_EQ(reference_param.surface_link().index(), 0u)
      << " Initial surface not found " << std::endl
      << " log10(RK tolerance): " << math::log10(rk_tolerance)
      << " Path length [mm]: " << bound_getter.m_path_length
      << " Average step size [mm]: " << bound_getter.m_avg_step_size
      << " Phi: " << reference_param.phi()
      << " Theta: " << reference_param.theta()
      << " Mom [GeV/c]: " << reference_param.p(ptc.charge());
  ASSERT_EQ(final_param.surface_link().index(), 1u)
      << " Final surface not found " << std::endl
      << " log10(RK tolerance): " << math::log10(rk_tolerance)
      << " Path length [mm]: " << bound_getter.m_path_length
      << " Average step size [mm]: " << bound_getter.m_avg_step_size
      << " Phi: " << reference_param.phi()
      << " Theta: " << reference_param.theta()
      << " Mom [GeV/c]: " << reference_param.p(ptc.charge());
  ASSERT_TRUE(detector_length > 0.f);
  ASSERT_GE(bound_getter.m_path_length, 0.5f * detector_length);
  ASSERT_LE(bound_getter.m_path_length, 1.5f * detector_length);
  ASSERT_LE(bound_getter.m_path_length, max_detector_length + 200.f);
  ASSERT_GE(bound_getter.m_abs_path_length, bound_getter.m_path_length);

  const auto reference_jacobian = bound_getter.m_jacobi;

  file << trk_count << ",";

  file << euler_angles_I[0u] << "," << euler_angles_I[1u] << ","
       << euler_angles_I[2u] << ",";
  file << euler_angles_F[0u] << "," << euler_angles_F[1u] << ","
       << euler_angles_F[2u] << ",";

  file << reference_param.bound_local()[0] << ","
       << reference_param.bound_local()[1] << "," << reference_param.phi()
       << "," << reference_param.theta() << "," << reference_param.qop() << ",";

  file << final_param.bound_local()[0] << "," << final_param.bound_local()[1]
       << "," << final_param.phi() << "," << final_param.theta() << ","
       << final_param.qop() << ",";

  bound_covariance_type differentiated_jacobian;
  std::array<unsigned int, 5u> num_iterations;
  std::array<bool, 25u> convergence;

  if (use_precal_values) {
    differentiated_jacobian = precal_diff_jacobi;
    num_iterations = precal_num_iterations;
    convergence = precal_convergence;
  } else {
    differentiated_jacobian = directly_differentiate<propagator_t, field_t>(
        trk_count, reference_param, det, detector_length, field,
        overstep_tolerance, path_tolerance, rk_tolerance_dis, constraint_step,
        hs, num_iterations, convergence);
  }

  bool total_convergence = (std::ranges::count(convergence, false) == 0);

  // Ridders number of iterations
  for (unsigned int i = 0; i < 5u; i++) {
    file << num_iterations[i] << ",";
  }

  // Convergence
  file << total_convergence << ",";
  for (unsigned int i = 0; i < 25u; i++) {
    file << convergence[i] << ",";
  }

  // Reference track
  for (unsigned int i = 0; i < 5u; i++) {
    for (unsigned int j = 0; j < 5u; j++) {
      file << getter::element(reference_jacobian, i, j) << ",";
    }
  }

  // Numerical evaluation
  for (unsigned int i = 0; i < 5u; i++) {
    for (unsigned int j = 0; j < 5u; j++) {
      file << getter::element(differentiated_jacobian, i, j) << ",";
    }
  }

  // Difference between evaluation and direct jacobian
  for (unsigned int i = 0; i < 5u; i++) {
    for (unsigned int j = 0; j < 5u; j++) {
      const scalar ref_val = getter::element(reference_jacobian, i, j);
      const scalar num_val = getter::element(differentiated_jacobian, i, j);
      const scalar rel_diff = get_relative_difference(ref_val, num_val);

      file << rel_diff << ",";

      // We return dqopdqop for test
      if (i == 4 && j == 4) {
        ref_rel_diff = rel_diff;
      }
    }
  }

  // Path length
  file << bound_getter.m_path_length << ",";

  // Absolute path length
  file << bound_getter.m_abs_path_length << ",";

  // The number of steps made
  file << bound_getter.step_count << ",";

  // Average step size
  file << bound_getter.m_avg_step_size << ",";

  // Log10(RK tolerance)
  file << math::log10(rk_tolerance) << ",";

  // Log10(on surface tolerance)
  file << math::log10(path_tolerance) << ",";

  // Overstep tolerance
  file << overstep_tolerance << ",";

  // Detector length
  file << detector_length;

  file << std::endl;
}

template <typename propagator_t, typename field_t>
void evaluate_covariance_transport(
    const std::size_t trk_count, const std::array<scalar, 3u>& euler_angles_I,
    const std::array<scalar, 3u>& euler_angles_F,
    const typename propagator_t::detector_type& det,
    const scalar detector_length,
    const bound_track_parameters<test_algebra>& track, const field_t& field,
    const scalar overstep_tolerance, const scalar path_tolerance,
    const scalar rk_tolerance, const scalar rk_tolerance_dis,
    const scalar constraint_step, std::ofstream& file,
    bool use_field_gradient) {
  // Copy track
  auto track_copy = track;

  // Make initial covariance
  const bound_covariance_type ini_cov =
      get_random_initial_covariance(track_copy.qop());

  track_copy.set_covariance(ini_cov);

  auto bound_getter = evaluate_bound_param<propagator_t, field_t>(
      trk_count, detector_length, track_copy, det, field, overstep_tolerance,
      path_tolerance, rk_tolerance, constraint_step, use_field_gradient, true,
      false);

  const auto reference_param = bound_getter.m_param_departure;
  const auto final_param = bound_getter.m_param_destination;
  const auto fin_cov = final_param.covariance();

  // Sanity check
  ASSERT_EQ(reference_param.surface_link().index(), 0u)
      << " Initial surface not found " << std::endl
      << " log10(RK tolerance): " << math::log10(rk_tolerance)
      << " Path length [mm]: " << bound_getter.m_path_length
      << " Average step size [mm]: " << bound_getter.m_avg_step_size
      << " Phi: " << reference_param.phi()
      << " Theta: " << reference_param.theta()
      << " Mom [GeV/c]: " << reference_param.p(ptc.charge());
  ASSERT_EQ(final_param.surface_link().index(), 1u)
      << " Final surface not found " << std::endl
      << " log10(RK tolerance): " << math::log10(rk_tolerance)
      << " Path length [mm]: " << bound_getter.m_path_length
      << " Average step size [mm]: " << bound_getter.m_avg_step_size
      << " Phi: " << reference_param.phi()
      << " Theta: " << reference_param.theta()
      << " Mom [GeV/c]: " << reference_param.p(ptc.charge());
  ASSERT_TRUE(detector_length > 0.f);
  ASSERT_GE(bound_getter.m_path_length, 0.5f * detector_length);
  ASSERT_LE(bound_getter.m_path_length, 1.5f * detector_length);
  ASSERT_LE(bound_getter.m_path_length, max_detector_length + 200.f);
  ASSERT_GE(bound_getter.m_abs_path_length, bound_getter.m_path_length);

  // Get smeared initial bound vector
  const bound_param_vector_type smeared_ini_vec =
      get_smeared_bound_vector(ini_cov, reference_param);

  // Make smeared bound track parameter
  auto smeared_track = track_copy;
  smeared_track.set_parameter_vector(smeared_ini_vec);

  auto smeared_bound_getter = evaluate_bound_param<propagator_t, field_t>(
      trk_count, detector_length, smeared_track, det, field, overstep_tolerance,
      path_tolerance, rk_tolerance_dis, constraint_step, use_field_gradient,
      false, false);

  // Get smeared final bound vector
  bound_param_vector_type smeared_fin_vec =
      smeared_bound_getter.m_param_destination;

  // phi needs to be wrapped w.r.t. phi of the reference vector
  wrap_angles(final_param, smeared_fin_vec);

  // Get pull values
  std::array<scalar, 5u> pulls;

  bound_vector<test_algebra> diff{};
  for (unsigned int i = 0u; i < 5u; i++) {
    getter::element(diff, i, 0u) = smeared_fin_vec[i] - final_param[i];
    pulls[i] = getter::element(diff, i, 0u) /
               math::sqrt(getter::element(fin_cov, i, i));
  }

  // Get Chi2
  const matrix_type<1u, 1u> chi2 =
      matrix::transpose(diff) * matrix::inverse(fin_cov) * diff;
  const scalar chi2_val = getter::element(chi2, 0u, 0u);

  file << trk_count << ",";

  file << euler_angles_I[0u] << "," << euler_angles_I[1u] << ","
       << euler_angles_I[2u] << ",";
  file << euler_angles_F[0u] << "," << euler_angles_F[1u] << ","
       << euler_angles_F[2u] << ",";

  // File writing
  file << reference_param[e_bound_loc0] << "," << reference_param[e_bound_loc1]
       << "," << reference_param[e_bound_phi] << ","
       << reference_param[e_bound_theta] << ","
       << reference_param[e_bound_qoverp] << ",";

  for (unsigned int i = 0; i < 5u; i++) {
    for (unsigned int j = 0; j < 5u; j++) {
      file << getter::element(ini_cov, i, j) << ",";
    }
  }

  file << final_param[e_bound_loc0] << "," << final_param[e_bound_loc1] << ","
       << final_param[e_bound_phi] << "," << final_param[e_bound_theta] << ","
       << final_param[e_bound_qoverp] << ",";

  for (unsigned int i = 0; i < 5u; i++) {
    for (unsigned int j = 0; j < 5u; j++) {
      file << getter::element(fin_cov, i, j) << ",";
    }
  }

  file << smeared_ini_vec[e_bound_loc0] << "," << smeared_ini_vec[e_bound_loc1]
       << "," << smeared_ini_vec[e_bound_phi] << ","
       << smeared_ini_vec[e_bound_theta] << ","
       << smeared_ini_vec[e_bound_qoverp] << ",";

  file << smeared_fin_vec[e_bound_loc0] << "," << smeared_fin_vec[e_bound_loc1]
       << "," << smeared_fin_vec[e_bound_phi] << ","
       << smeared_fin_vec[e_bound_theta] << ","
       << smeared_fin_vec[e_bound_qoverp] << ",";

  file << pulls[0] << "," << pulls[1] << "," << pulls[2] << "," << pulls[3]
       << "," << pulls[4] << ",";

  file << chi2_val << ",";

  // Path length
  file << bound_getter.m_path_length << ",";

  // Absolute path length
  file << bound_getter.m_abs_path_length << ",";

  // The number of steps made
  file << bound_getter.step_count << ",";

  // Average step size
  file << bound_getter.m_avg_step_size << ",";

  // Log10(RK tolerance)
  file << math::log10(rk_tolerance) << ",";

  // Log10(on surface tolerance)
  file << math::log10(path_tolerance) << ",";

  // Overstep tolerance
  file << overstep_tolerance << ",";

  // Detector length
  file << detector_length;

  file << std::endl;
}

template <typename detector_t, typename detector_t::metadata::mask_id mask_id>
bound_param_vector_type get_displaced_bound_vector_helix(
    const bound_track_parameters<test_algebra>& track, const vector3& field,
    unsigned int target_index, scalar displacement, const detector_t& det,
    const scalar helix_tolerance) {
  const auto& departure_sf = det.surface(0u);

  const auto& destination_sf = det.surface(1u);
  const auto& trf_link = destination_sf.transform();
  const auto& destination_trf = det.transform_store().at(trf_link);
  const auto& mask_link = destination_sf.mask();
  const auto& destination_mask =
      det.mask_store().template get<mask_id>().at(mask_link.index());

  bound_param_vector_type dvec = track;
  dvec[target_index] += displacement;
  const auto free_vec =
      tracking_surface{det, departure_sf}.bound_to_free_vector({}, dvec);
  detail::helix<test_algebra> hlx(free_vec, field);

  using mask_t = types::get<typename detector_t::masks, mask_id>;
  helix_intersector<typename mask_t::shape, test_algebra> hlx_is{};
  hlx_is.run_rtsafe = false;
  hlx_is.convergence_tolerance = helix_tolerance;
  auto sfi =
      hlx_is(hlx, destination_sf, destination_mask, destination_trf, 0.f);
  const auto path_length = sfi.path();
  const auto pos = hlx(path_length);
  const auto dir = hlx.dir(path_length);

  const free_track_parameters<test_algebra> new_free_par(pos, 0, dir,
                                                         hlx.qop());
  auto new_bound_vec =
      tracking_surface{det, destination_sf}.free_to_bound_vector({},
                                                                 new_free_par);

  // phi needs to be wrapped w.r.t. phi of the reference vector
  wrap_angles(dvec, new_bound_vec);

  return new_bound_vec;
}

template <typename detector_t, typename detector_t::metadata::mask_id mask_id>
void evaluate_jacobian_difference_helix(
    const std::size_t trk_count, const std::array<scalar, 3u> euler_angles_I,
    const std::array<scalar, 3u> euler_angles_F, const detector_t& det,
    const scalar detector_length,
    const bound_track_parameters<test_algebra>& track, const vector3& field,
    const std::array<scalar, 5u> hs, std::ofstream& file,
    const scalar helix_tolerance) {
  const auto phi0 = track.phi();
  const auto theta0 = track.theta();
  (void)phi0;
  (void)theta0;

  // Get bound to free Jacobi
  const auto& departure_sf = det.surface(0u);
  const auto bound_to_free_jacobi =
      tracking_surface{det, departure_sf}.bound_to_free_jacobian({}, track);

  // Get fre vector
  const auto free_vec =
      tracking_surface{det, departure_sf}.bound_to_free_vector({}, track);

  // Helix from the departure surface
  detail::helix<test_algebra> hlx(free_vec, field);

  const auto& destination_sf = det.surface(1u);
  const auto& trf_link = destination_sf.transform();
  const auto& destination_trf = det.transform_store().at(trf_link);
  const auto& mask_link = destination_sf.mask();
  const auto& destination_mask =
      det.mask_store().template get<mask_id>().at(mask_link.index());

  using mask_t = types::get<typename detector_t::masks, mask_id>;
  helix_intersector<typename mask_t::shape, test_algebra> hlx_is{};
  hlx_is.run_rtsafe = false;
  hlx_is.convergence_tolerance = helix_tolerance;

  auto sfi =
      hlx_is(hlx, destination_sf, destination_mask, destination_trf, 0.f);

  EXPECT_TRUE(sfi.is_inside())
      << " Final surface not found" << std::endl
      << " log10(Helix tolerance): " << math::log10(helix_tolerance)
      << " Phi: " << track.phi() << " Theta: " << track.theta()
      << " Mom [GeV/c]: " << track.p(ptc.charge());

  const auto path_length = sfi.path();

  // Get transport Jacobi
  const auto transport_jacobi = hlx.jacobian(path_length);

  const auto pos = hlx(path_length);
  const auto dir = hlx.dir(path_length);
  const auto qop = hlx.qop();

  // Get correction term
  const auto correction_term =
      matrix::identity<free_matrix<test_algebra>>() +
      tracking_surface{det, destination_sf}.path_correction(
          {}, pos, dir, qop * vector::cross(dir, field), 0.f);

  const free_track_parameters<test_algebra> free_par(pos, 0.f, dir, qop);

  // Get free to bound Jacobi
  const auto free_to_bound_jacobi =
      tracking_surface{det, destination_sf}.free_to_bound_jacobian({},
                                                                   free_par);

  // Get full Jacobi
  const auto reference_jacobian = free_to_bound_jacobi * correction_term *
                                  transport_jacobi * bound_to_free_jacobi;

  // Get bound vector
  const auto bound_vec =
      tracking_surface{det, destination_sf}.free_to_bound_vector({}, free_par);

  /******************************
   *  Numerical differentiation
   * ****************************/

  bound_covariance_type differentiated_jacobian;
  std::array<unsigned int, 5u> num_iterations;
  std::array<bool, 25u> convergence;

  for (unsigned int i = 0u; i < 5u; i++) {
    scalar delta = hs[i];

    preprocess_delta(i, delta, track);

    const auto vec1 = get_displaced_bound_vector_helix<detector_t, mask_id>(
        track, field, i, 1.f * delta, det, helix_tolerance);
    const auto vec2 = get_displaced_bound_vector_helix<detector_t, mask_id>(
        track, field, i, -1.f * delta, det, helix_tolerance);

    ridders_derivative ridder;
    ridder.initialize(vec1, vec2, delta);

    for (unsigned int p = 1u; p < Nt; p++) {
      delta /= con[i];

      const auto nvec1 = get_displaced_bound_vector_helix<detector_t, mask_id>(
          track, field, i, 1.f * delta, det, helix_tolerance);
      const auto nvec2 = get_displaced_bound_vector_helix<detector_t, mask_id>(
          track, field, i, -1.f * delta, det, helix_tolerance);

      ridder.run(nvec1, nvec2, delta, p, i, differentiated_jacobian);

      if (ridder.finished()) {
        num_iterations[i] = p;
        break;
      }
    }

    for (std::size_t j = 0u; j < 5u; j++) {
      convergence[i + j * 5] = ridder.complete[j];
    }
  }

  bool total_convergence = (std::ranges::count(convergence, false) == 0);

  file << trk_count << ",";

  file << euler_angles_I[0u] << "," << euler_angles_I[1u] << ","
       << euler_angles_I[2u] << ",";
  file << euler_angles_F[0u] << "," << euler_angles_F[1u] << ","
       << euler_angles_F[2u] << ",";

  file << track.bound_local()[0] << "," << track.bound_local()[1] << ","
       << track.phi() << "," << track.theta() << "," << track.qop() << ",";

  file << bound_vec[e_bound_loc0] << "," << bound_vec[e_bound_loc1] << ","
       << bound_vec.phi() << "," << bound_vec.theta() << "," << bound_vec.qop()
       << ",";

  // Ridders number of iterations
  for (unsigned int i = 0; i < 5u; i++) {
    file << num_iterations[i] << ",";
  }

  // Convergence
  file << total_convergence << ",";
  for (unsigned int i = 0; i < 25u; i++) {
    file << convergence[i] << ",";
  }

  // Reference track
  for (unsigned int i = 0; i < 5u; i++) {
    for (unsigned int j = 0; j < 5u; j++) {
      file << getter::element(reference_jacobian, i, j) << ",";
    }
  }

  // Numerical evaluation
  for (unsigned int i = 0; i < 5u; i++) {
    for (unsigned int j = 0; j < 5u; j++) {
      file << getter::element(differentiated_jacobian, i, j) << ",";
    }
  }

  // Difference between evaluation and direct jacobian
  for (unsigned int i = 0; i < 5u; i++) {
    for (unsigned int j = 0; j < 5u; j++) {
      const scalar ref_val = getter::element(reference_jacobian, i, j);
      const scalar num_val = getter::element(differentiated_jacobian, i, j);
      const scalar rel_diff = get_relative_difference(ref_val, num_val);
      file << rel_diff << ",";
    }
  }

  // Path length
  file << path_length << ",";

  // Absolute Path length (Not supported for helix)
  file << 0 << ",";

  // Average step size (Doesn't exist for helix intersection)
  file << 0 << ",";

  // The number of steps made (Doesn't exist for helix intersection)
  file << 0 << ",";

  // Log10(RK tolerance) (Doesn't exist for helix intersection)
  file << 0 << ",";

  // Log10(helix intersection tolerance)
  file << math::log10(helix_tolerance) << ",";

  // Overstep tolerance (Doesn't exist for helix intersection)
  file << 0 << ",";

  // Detector length
  file << detector_length;

  file << std::endl;
}

void setup_csv_header_jacobian(std::ofstream& file) {
  file << std::fixed << std::showpoint;
  file << std::setprecision(32);

  // Track ID
  file << "track_ID,";

  // Euler angles
  file << "alpha_I,beta_I,gamma_I,";
  file << "alpha_F,beta_F,gamma_F,";

  // Initial Parameter at the departure surface
  file << "l0_I,l1_I,phi_I,theta_I,qop_I,";

  // Final Parameter at the destination surface
  file << "l0_F,l1_F,phi_F,theta_F,qop_F,";

  // Number of iterations to complete the numerical differentiation
  file << "num_iterations_l0,"
       << "num_iterations_l1,"
       << "num_iterations_phi,"
       << "num_iterations_theta,"
       << "num_iterations_qop,";

  // Convergence
  file << "total_convergence,";
  file << "dl0dl0_C,dl0dl1_C,dl0dphi_C,dl0dtheta_C,dl0dqop_C,";
  file << "dl1dl0_C,dl1dl1_C,dl1dphi_C,dl1dtheta_C,dl1dqop_C,";
  file << "dphidl0_C,dphidl1_C,dphidphi_C,dphidtheta_C,dphidqop_C,";
  file << "dthetadl0_C,dthetadl1_C,dthetadphi_C,dthetadtheta_C,"
          "dthetadqop_C,";
  file << "dqopdl0_C,dqopdl1_C,dqopdphi_C,dqopdtheta_C,dqopdqop_C,";

  // Evaluation
  file << "dl0dl0_E,dl0dl1_E,dl0dphi_E,dl0dtheta_E,dl0dqop_E,";
  file << "dl1dl0_E,dl1dl1_E,dl1dphi_E,dl1dtheta_E,dl1dqop_E,";
  file << "dphidl0_E,dphidl1_E,dphidphi_E,dphidtheta_E,dphidqop_E,";
  file << "dthetadl0_E,dthetadl1_E,dthetadphi_E,dthetadtheta_E,"
          "dthetadqop_E,";
  file << "dqopdl0_E,dqopdl1_E,dqopdphi_E,dqopdtheta_E,dqopdqop_E,";

  // Numerical Differentiation
  file << "dl0dl0_D,dl0dl1_D,dl0dphi_D,dl0dtheta_D,dl0dqop_D,";
  file << "dl1dl0_D,dl1dl1_D,dl1dphi_D,dl1dtheta_D,dl1dqop_D,";
  file << "dphidl0_D,dphidl1_D,dphidphi_D,dphidtheta_D,dphidqop_D,";
  file << "dthetadl0_D,dthetadl1_D,dthetadphi_D,dthetadtheta_D,"
          "dthetadqop_D,";
  file << "dqopdl0_D,dqopdl1_D,dqopdphi_D,dqopdtheta_D,dqopdqop_D,";

  // Relative Difference between Evaluation and Numerical Differentiation
  file << "dl0dl0_R,dl0dl1_R,dl0dphi_R,dl0dtheta_R,dl0dqop_R,";
  file << "dl1dl0_R,dl1dl1_R,dl1dphi_R,dl1dtheta_R,dl1dqop_R,";
  file << "dphidl0_R,dphidl1_R,dphidphi_R,dphidtheta_R,dphidqop_R,";
  file << "dthetadl0_R,dthetadl1_R,dthetadphi_R,dthetadtheta_R,"
          "dthetadqop_R,";
  file << "dqopdl0_R,dqopdl1_R,dqopdphi_R,dqopdtheta_R,dqopdqop_R,";

  // Path length [mm]
  file << "path_length,";

  // Absolute path length [mm]
  file << "abs_path_length,";

  // The number of steps made
  file << "n_steps,";

  // Average step size [mm]
  file << "average_step_size,";

  // RK Tolerances [mm]
  file << "log10_rk_tolerance,";

  // Intersection tolerance [mm]
  file << "log10_intersection_tolerance,";

  // Overstep tolerance [mm]
  file << "overstep_tolerance,";

  // Detector length [mm]
  file << "detector_length";

  file << std::endl;
}

void setup_csv_header_covariance(std::ofstream& file) {
  file << std::fixed << std::showpoint;
  file << std::setprecision(32);

  // Track ID
  file << "track_ID,";

  // Euler angles
  file << "alpha_I,beta_I,gamma_I,";
  file << "alpha_F,beta_F,gamma_F,";

  // Initial parameters (vector + covariance) of reference track
  file << "l0_I,l1_I,phi_I,theta_I,qop_I,";
  file << "l0l0_I,l0l1_I,l0phi_I,l0theta_I,l0qop_I,";
  file << "l1l0_I,l1l1_I,l1phi_I,l1theta_I,l1qop_I,";
  file << "phil0_I,phil1_I,phiphi_I,phitheta_I,phiqop_I,";
  file << "thetal0_I,thetal1_I,thetaphi_I,thetatheta_I,thetaqop_I,";
  file << "qopl0_I,qopl1_I,qopphi_I,qoptheta_I,qopqop_I,";

  // Final parameters (vector + covariance) of reference track
  file << "l0_F,l1_F,phi_F,theta_F,qop_F,";
  file << "l0l0_F,l0l1_F,l0phi_F,l0theta_F,l0qop_F,";
  file << "l1l0_F,l1l1_F,l1phi_F,l1theta_F,l1qop_F,";
  file << "phil0_F,phil1_F,phiphi_F,phitheta_F,phiqop_F,";
  file << "thetal0_F,thetal1_F,thetaphi_F,thetatheta_F,thetaqop_F,";
  file << "qopl0_F,qopl1_F,qopphi_F,qoptheta_F,qopqop_F,";

  // Initial parameter (vector only) of smeared track
  file << "l0_IS,l1_IS,phi_IS,theta_IS,qop_IS,";

  // Final parameter (vector only) of smeared track
  file << "l0_FS,l1_FS,phi_FS,theta_FS,qop_FS,";

  // Pull values
  file << "pull_l0,pull_l1,pull_phi,pull_theta,pull_qop,";

  // Chi2
  file << "chi2,";

  // Path length [mm]
  file << "path_length,";

  // Absolute path length [mm]
  file << "abs_path_length,";

  // The number of steps made
  file << "n_steps,";

  // Average step size [mm]
  file << "average_step_size,";

  // Tolerances [mm]
  file << "log10_rk_tolerance,";

  // Intersection tolerance [mm]
  file << "log10_intersection_tolerance,";

  // Overstep tolerance [mm]
  file << "overstep_tolerance,";

  // Detector length [mm]
  file << "detector_length";

  file << std::endl;
}

int main(int argc, char** argv) {
  // Options parsing
  po::options_description desc("\ndetray jacobian validation options");
  desc.add_options()("help", "produce help message");
  desc.add_options()("output-directory",
                     po::value<std::string>()->default_value(""),
                     "Output directory");
  desc.add_options()("n-tracks", po::value<std::size_t>()->default_value(100u),
                     "Number of tracks for generator");
  desc.add_options()("n-skips", po::value<std::size_t>()->default_value(0u),
                     "Number of skipped indices");
  desc.add_options()("skip-rect", po::value<bool>()->default_value(false),
                     "Skip rectangular telescope");
  desc.add_options()("skip-wire", po::value<bool>()->default_value(false),
                     "Skip wire telescope");
  desc.add_options()("log10-rk-tolerance-jac-mm",
                     po::value<scalar>()->default_value(-4.f),
                     "Set log10(rk_tolerance_jac_in_mm)");
  desc.add_options()("log10-rk-tolerance-dis-mm",
                     po::value<scalar>()->default_value(-6.f),
                     "Set log10(rk_tolerance_dis_in_mm)");
  desc.add_options()("log10-rk-tolerance-cov-mm",
                     po::value<scalar>()->default_value(-4.f),
                     "Set log10(rk_tolerance_cov_in_mm)");
  desc.add_options()("log10-helix-tolerance-mm",
                     po::value<scalar>()->default_value(-3.f),
                     "Set log10(helix_tolerance_in_mm)");
  desc.add_options()("overstep-tolerance-mm",
                     po::value<scalar>()->default_value(-1000.f),
                     "Set the overstep tolerance in mm unit");
  desc.add_options()("log10-on-surface-tolerance-mm",
                     po::value<scalar>()->default_value(-3.f),
                     "Set log10(path_tolerance_in_mm)");
  desc.add_options()("rk-tolerance-iterate-mode",
                     po::value<bool>()->default_value(true),
                     "Iterate over the rk tolerances");
  desc.add_options()("log10-min-rk-tolerance-mm",
                     po::value<scalar>()->default_value(-6.f),
                     "Set log10(min_rk_tolerance_in_mm)");
  desc.add_options()("log10-max-rk-tolerance-mm",
                     po::value<scalar>()->default_value(2.f),
                     "Set log10(max_rk_tolerance_in_mm)");
  desc.add_options()("mc-seed", po::value<std::size_t>()->default_value(0u),
                     "Monte-Carlo seed");
  desc.add_options()("verbose-level", po::value<int>()->default_value(1),
                     "Verbose level");

  po::variables_map vm;
  po::store(parse_command_line(argc, argv, desc,
                               po::command_line_style::unix_style ^
                                   po::command_line_style::allow_short),
            vm);
  po::notify(vm);

  // Help message
  if (vm.count("help")) {
    std::clog << desc << std::endl;
    return EXIT_FAILURE;
  }

  const std::string output_directory = vm["output-directory"].as<std::string>();
  std::size_t n_tracks = vm["n-tracks"].as<std::size_t>();
  std::size_t n_skips = vm["n-skips"].as<std::size_t>();
  const bool skip_rect = vm["skip-rect"].as<bool>();
  const bool skip_wire = vm["skip-wire"].as<bool>();
  const scalar rk_power_jac = vm["log10-rk-tolerance-jac-mm"].as<scalar>();
  const scalar rk_tol_jac = std::pow(10.f, rk_power_jac) * unit<scalar>::mm;
  const scalar rk_power_dis = vm["log10-rk-tolerance-dis-mm"].as<scalar>();
  const scalar rk_tol_dis = std::pow(10.f, rk_power_dis) * unit<scalar>::mm;
  const scalar rk_power_cov = vm["log10-rk-tolerance-cov-mm"].as<scalar>();
  const scalar rk_tol_cov = std::pow(10.f, rk_power_cov) * unit<scalar>::mm;
  const scalar helix_power = vm["log10-helix-tolerance-mm"].as<scalar>();
  const scalar helix_tol = std::pow(10.f, helix_power) * unit<scalar>::mm;
  const scalar on_surface_power =
      vm["log10-on-surface-tolerance-mm"].as<scalar>() * unit<scalar>::mm;
  const scalar on_surface_tol = math::pow(10.f, on_surface_power);
  const bool rk_tolerance_iterate_mode =
      vm["rk-tolerance-iterate-mode"].as<bool>();
  const scalar log10_min_rk_tolerance =
      vm["log10-min-rk-tolerance-mm"].as<scalar>() * unit<scalar>::mm;
  const scalar log10_max_rk_tolerance =
      vm["log10-max-rk-tolerance-mm"].as<scalar>() * unit<scalar>::mm;
  const std::size_t mc_seed = vm["mc-seed"].as<std::size_t>();
  const int verbose_lvl = vm["verbose-level"].as<int>();

  std::vector<scalar> log10_tols;
  scalar r = log10_min_rk_tolerance;
  while (r <= log10_max_rk_tolerance + 1e-3f) {
    log10_tols.push_back(r);
    r = r + 2.f;
  }

  // Set seed for random generator
  mt1.seed(mc_seed);

  // Volume material
  const material<scalar> volume_mat = detray::cesium_iodide_with_ded<scalar>();

  std::string path;
  // Create output directory
  if (output_directory.empty()) {
    path = "";
  } else {
    std::filesystem::create_directories(output_directory);
    path = output_directory + "/";
  }

  // Output Csv file
  std::ofstream helix_rect_file;
  std::ofstream const_rect_file;
  std::ofstream inhom_rect_file;
  std::ofstream inhom_rect_material_file;
  std::ofstream helix_wire_file;
  std::ofstream const_wire_file;
  std::ofstream inhom_wire_file;
  std::ofstream inhom_wire_material_file;
  helix_rect_file.open(path + "helix_rect.csv");
  const_rect_file.open(path + "const_rect.csv");
  inhom_rect_file.open(path + "inhom_rect.csv");
  inhom_rect_material_file.open(path + "inhom_rect_material.csv");
  helix_wire_file.open(path + "helix_wire.csv");
  const_wire_file.open(path + "const_wire.csv");
  inhom_wire_file.open(path + "inhom_wire.csv");
  inhom_wire_material_file.open(path + "inhom_wire_material.csv");

  setup_csv_header_jacobian(helix_rect_file);
  setup_csv_header_jacobian(const_rect_file);
  setup_csv_header_jacobian(inhom_rect_file);
  setup_csv_header_jacobian(inhom_rect_material_file);
  setup_csv_header_jacobian(helix_wire_file);
  setup_csv_header_jacobian(const_wire_file);
  setup_csv_header_jacobian(inhom_wire_file);
  setup_csv_header_jacobian(inhom_wire_material_file);

  std::ofstream rect_cov_transport_file;
  rect_cov_transport_file.open(path + "rect_cov_transport.csv");
  setup_csv_header_covariance(rect_cov_transport_file);

  std::ofstream wire_cov_transport_file;
  wire_cov_transport_file.open(path + "wire_cov_transport.csv");
  setup_csv_header_covariance(wire_cov_transport_file);

  // Output Csv file (RK tolerance iteration mode)
  std::vector<std::ofstream> rect_files(log10_tols.size());
  std::vector<std::ofstream> wire_files(log10_tols.size());

  if (rk_tolerance_iterate_mode) {
    for (std::size_t i = 0u; i < log10_tols.size(); i++) {
      const std::string rect_name =
          "inhom_rect_material_" + std::to_string(int(log10_tols[i])) + ".csv";
      const std::string wire_name =
          "inhom_wire_material_" + std::to_string(int(log10_tols[i])) + ".csv";
      rect_files[i].open(path + rect_name);
      wire_files[i].open(path + wire_name);

      setup_csv_header_jacobian(rect_files[i]);
      setup_csv_header_jacobian(wire_files[i]);
    }
  }

  // Memory resource
  vecmem::host_memory_resource host_mr;

  // Filter out the google test flags
  ::testing::InitGoogleTest(&argc, argv);

  // Detector types
  using rectangle_telescope =
      detector<telescope_metadata<test_algebra, rect_type>>;
  using wire_telescope = detector<telescope_metadata<test_algebra, wire_type>>;
  using track_type = free_track_parameters<test_algebra>;

  // Constant magnetic field type
  using const_bfield_t = bfield::const_field_t<scalar>;

  // Magnetic field map using nearest neighbor interpolation
  using inhom_bfield_t = bfield::inhom_field_t<scalar>;

  const const_bfield_t const_bfield = create_const_field<scalar>(B_z);
  const inhom_bfield_t inhom_bfield = create_inhom_field<scalar>();

  // Actor chain type
  using actor_chain_t = actor_chain<
      actor::parameter_updater<test_algebra, bound_getter<test_algebra>>>;

  // Iterate over reference (pilot) tracks for a rectangular telescope
  // geometry and Jacobian calculation
  using uniform_gen_t =
      detail::random_numbers<scalar, std::uniform_real_distribution<scalar>>;
  using trk_generator_t = random_track_generator<track_type, uniform_gen_t>;
  trk_generator_t::configuration trk_gen_cfg{};
  trk_gen_cfg.seed(mc_seed);
  trk_gen_cfg.n_tracks(n_tracks + n_skips);
  trk_gen_cfg.phi_range(-constant<scalar>::pi, constant<scalar>::pi);
  trk_gen_cfg.theta_range(0.f, constant<scalar>::pi);
  trk_gen_cfg.mom_range(min_mom, max_mom);
  trk_gen_cfg.origin(0.f, 0.f, 0.f);
  trk_gen_cfg.origin_stddev(0.f, 0.f, 0.f);

  // Vectors for dqopdqop relative difference
  std::vector<std::vector<scalar>> dqopdqop_rel_diffs_rect(log10_tols.size());
  std::vector<std::vector<scalar>> dqopdqop_rel_diffs_wire(log10_tols.size());

  bool do_inspect = false;
  if (verbose_lvl >= 4) {
    do_inspect = true;
  }

  // Navigator types
  using rect_navigator_t = caching_navigator<rectangle_telescope>;
  using wire_navigator_t = caching_navigator<wire_telescope>;

  // Stepper types
  using const_field_stepper_t =
      rk_stepper<const_bfield_t::view_t, test_algebra, constrained_step<scalar>,
                 stepper_default_policy<scalar>>;
  using inhom_field_stepper_t =
      rk_stepper<inhom_bfield_t::view_t, test_algebra, constrained_step<scalar>,
                 stepper_default_policy<scalar>>;

  // Make four propagators for each case
  using const_field_rect_propagator_t =
      propagator<const_field_stepper_t, rect_navigator_t, actor_chain_t>;
  using inhom_field_rect_propagator_t =
      propagator<inhom_field_stepper_t, rect_navigator_t, actor_chain_t>;

  using const_field_wire_propagator_t =
      propagator<const_field_stepper_t, wire_navigator_t, actor_chain_t>;
  using inhom_field_wire_propagator_t =
      propagator<inhom_field_stepper_t, wire_navigator_t, actor_chain_t>;

  std::size_t track_count = 0u;

  for (const auto track : trk_generator_t{trk_gen_cfg}) {
    mt2.seed(track_count);

    // Pilot track
    detail::helix<test_algebra> helix_bz(track, B_z);

    // Make a telescope geometry with rectagular surface
    const scalar detector_length = rand_length(mt1);
    const scalar constraint_step_size = detector_length * 1.25f;

    mask<rect_type, test_algebra> rect{0u, detector_length * mask_scaler,
                                       detector_length * mask_scaler};
    mask<wire_type, test_algebra> wire{0u, detector_length * mask_scaler,
                                       detector_length * mask_scaler};

    // Adjust overstep tolerance
    scalar overstep_tol =
        vm["overstep-tolerance-mm"].as<scalar>() * unit<scalar>::mm;
    overstep_tol = -math::min(math::fabs(overstep_tol), detector_length * 0.5f);

    tel_det_config<test_algebra, rect_type, detail::helix> rectangle_cfg{
        rect, helix_bz};
    rectangle_cfg.envelope(envelope_size);
    rectangle_cfg.module_material(vacuum<scalar>{});
    rectangle_cfg.mat_thickness(0.f);
    rectangle_cfg.n_surfaces(2u);
    rectangle_cfg.length(detector_length);
    rectangle_cfg.volume_material(vacuum<scalar>{});
    rectangle_cfg.do_check(false);

    auto alphaI = rand_alpha(mt1);
    auto alphaF = rand_alpha(mt1);
    auto betaF = math::acos(rand_cosbeta(mt1));
    if (rand_bool(mt1) == 0) {
      betaF = -betaF;
    }
    auto gammaF = rand_gamma(mt1);

    std::array<scalar, 3u> euler_angles_I{alphaI, 0.f, 0.f};
    std::array<scalar, 3u> euler_angles_F{alphaF, betaF, gammaF};

    // Without volume material
    auto [rect_det, rect_names] =
        build_telescope_detector<test_algebra>(host_mr, rectangle_cfg);
    const auto [euler_rect_initial, shift_rect_initial] =
        tilt_surface<decltype(rect_det),
                     decltype(rect_det)::masks::id::e_rectangle2D>(
            rect_det, 0u, helix_bz.dir(0.f), alphaI, 0.f, 0.f);
    const auto [euler_rect_final, shift_rect_final] =
        tilt_surface<decltype(rect_det),
                     decltype(rect_det)::masks::id::e_rectangle2D>(
            rect_det, 1u, helix_bz.dir(detector_length), alphaF, betaF, gammaF);

    // With volume material
    rectangle_cfg.volume_material(volume_mat);
    auto [rect_det_w_mat, rect_names2] =
        build_telescope_detector<test_algebra>(host_mr, rectangle_cfg);
    [[maybe_unused]] const auto [euler_rect_initial2, shift_rect_initial2] =
        tilt_surface<decltype(rect_det_w_mat),
                     decltype(rect_det_w_mat)::masks::id::e_rectangle2D>(
            rect_det_w_mat, 0u, helix_bz.dir(0.f), alphaI, 0.f, 0.f);
    [[maybe_unused]] const auto [euler_rect_final2, shift_rect_final2] =
        tilt_surface<decltype(rect_det_w_mat),
                     decltype(rect_det_w_mat)::masks::id::e_rectangle2D>(
            rect_det_w_mat, 1u, helix_bz.dir(detector_length), alphaF, betaF,
            gammaF);

    // Make a telescope geometry with wire surface
    tel_det_config<test_algebra, wire_type, detail::helix> wire_cfg{wire,
                                                                    helix_bz};
    wire_cfg.envelope(envelope_size);
    wire_cfg.module_material(vacuum<scalar>{});
    wire_cfg.mat_thickness(0.f);
    wire_cfg.n_surfaces(2u);
    wire_cfg.length(detector_length);
    wire_cfg.volume_material(vacuum<scalar>{});
    wire_cfg.do_check(false);

    // Without volume material
    auto [wire_det, wire_names] =
        build_telescope_detector<test_algebra>(host_mr, wire_cfg);
    const auto [euler_wire_initial, shift_wire_initial] =
        tilt_surface<decltype(wire_det),
                     decltype(wire_det)::masks::id::e_drift_cell>(
            wire_det, 0u, helix_bz.dir(0.f), alphaI, 0.f, 0.f);
    const auto [euler_wire_final, shift_wire_final] =
        tilt_surface<decltype(wire_det),
                     decltype(wire_det)::masks::id::e_drift_cell>(
            wire_det, 1u, helix_bz.dir(detector_length), alphaF, betaF, gammaF);

    // With volume material
    wire_cfg.volume_material(volume_mat);
    auto [wire_det_w_mat, wire_names2] =
        build_telescope_detector<test_algebra>(host_mr, wire_cfg);
    [[maybe_unused]] const auto [euler_wire_initial2, shift_wire_initial2] =
        tilt_surface<decltype(wire_det_w_mat),
                     decltype(wire_det_w_mat)::masks::id::e_drift_cell>(
            wire_det_w_mat, 0u, helix_bz.dir(0.f), alphaI, 0.f, 0.f);
    [[maybe_unused]] const auto [euler_wire_final2, shift_wire_final2] =
        tilt_surface<decltype(wire_det_w_mat),
                     decltype(wire_det_w_mat)::masks::id::e_drift_cell>(
            wire_det_w_mat, 1u, helix_bz.dir(detector_length), alphaF, betaF,
            gammaF);

    // This IF block should locate after `tilt_surface()` calls for
    // debugging purpose
    if (track_count + 1 <= n_skips) {
      track_count++;
      continue;
    }

    track_count++;

    if (verbose_lvl >= 1) {
      std::clog << "[Event Property]" << std::endl;
      std::clog << "Track ID: " << track_count
                << "  Number of processed tracks per thread: "
                << track_count - n_skips << std::endl;
    }
    if (verbose_lvl >= 2) {
      std::clog << "[Detector Property]" << std::endl;
      std::clog << "Path length for the final surface: " << detector_length
                << std::endl;
      std::clog << "Rect initial surface rotation: ("
                << euler_rect_initial.alpha << " " << euler_rect_initial.beta
                << " " << euler_rect_initial.gamma << ")" << std::endl;
      std::clog << "Rect initial surface shift: (" << shift_rect_initial[0u]
                << " " << shift_rect_initial[1u] << " "
                << shift_rect_initial[2u] << ")" << std::endl;
      std::clog << "Rect final surface rotation: (" << euler_rect_final.alpha
                << " " << euler_rect_final.beta << " " << euler_rect_final.gamma
                << ")" << std::endl;
      std::clog << "Rect final surface shift: (" << shift_rect_final[0u] << " "
                << shift_rect_final[1u] << " " << shift_rect_final[2u] << ")"
                << std::endl;
      std::clog << "Wire initial surface rotation: ("
                << euler_wire_initial.alpha << " " << euler_wire_initial.beta
                << " " << euler_wire_initial.gamma << ")" << std::endl;
      std::clog << "Wire initial surface shift: (" << shift_wire_initial[0u]
                << " " << shift_wire_initial[1u] << " "
                << shift_wire_initial[2u] << ")" << std::endl;
      std::clog << "Wire final surface rotation: (" << euler_wire_final.alpha
                << " " << euler_wire_final.beta << " " << euler_wire_final.gamma
                << ")" << std::endl;
      std::clog << "Wire final surface shift: (" << shift_wire_final[0u] << " "
                << shift_wire_final[1u] << " " << shift_wire_final[2u] << ")"
                << std::endl;
    }
    if (verbose_lvl >= 3) {
      std::clog << "[Track Property]" << std::endl;
      std::clog << "Phi: " << vector::phi(track.dir()) << std::endl;
      std::clog << "Theta: " << vector::theta(track.dir()) << std::endl;
      std::clog << "Mom: " << track.p(ptc.charge()) << std::endl;
    }

    /**********************************
     * Rectangular telescope geometry
     **********************************/

    if (verbose_lvl >= 3 && !skip_rect) {
      std::clog << "Simulating rectangular telescope..." << std::endl;
    }

    // Get initial parameter
    const auto rect_bparam =
        get_initial_parameter<decltype(rect_det),
                              decltype(rect_det)::masks::id::e_rectangle2D>(
            rect_det, track, B_z, helix_tol);

    if (!skip_rect) {
      scalar ref_rel_diff;

      if (rk_tolerance_iterate_mode) {
        std::array<unsigned int, 5u> num_iterations;
        std::array<bool, 25u> convergence;
        auto differentiated_jacobian =
            directly_differentiate<inhom_field_rect_propagator_t>(
                track_count, rect_bparam, rect_det_w_mat, detector_length,
                inhom_bfield, overstep_tol, on_surface_tol, rk_tol_dis,
                constraint_step_size, h_sizes_rect, num_iterations,
                convergence);

        for (std::size_t i = 0u; i < log10_tols.size(); i++) {
          // Rectangle Inhomogeneous field with Material
          evaluate_jacobian_difference<inhom_field_rect_propagator_t>(
              track_count, euler_angles_I, euler_angles_F, rect_det_w_mat,
              detector_length, rect_bparam, inhom_bfield, overstep_tol,
              on_surface_tol, std::pow(10.f, log10_tols[i]), rk_tol_dis,
              constraint_step_size, h_sizes_rect, rect_files[i], ref_rel_diff,
              true, do_inspect, true, differentiated_jacobian, num_iterations,
              convergence);

          dqopdqop_rel_diffs_rect[i].push_back(ref_rel_diff);
        }
      } else if (!rk_tolerance_iterate_mode) {
        // For helix
        evaluate_jacobian_difference_helix<
            decltype(rect_det), decltype(rect_det)::masks::id::e_rectangle2D>(
            track_count, euler_angles_I, euler_angles_F, rect_det,
            detector_length, rect_bparam, B_z, h_sizes_rect, helix_rect_file,
            helix_tol);

        // Rect Const field
        evaluate_jacobian_difference<const_field_rect_propagator_t>(
            track_count, euler_angles_I, euler_angles_F, rect_det,
            detector_length, rect_bparam, const_bfield, overstep_tol,
            on_surface_tol, rk_tol_jac, rk_tol_dis, constraint_step_size,
            h_sizes_rect, const_rect_file, ref_rel_diff, true, false);

        // Rect Inhomogeneous field
        evaluate_jacobian_difference<inhom_field_rect_propagator_t>(
            track_count, euler_angles_I, euler_angles_F, rect_det,
            detector_length, rect_bparam, inhom_bfield, overstep_tol,
            on_surface_tol, rk_tol_jac, rk_tol_dis, constraint_step_size,
            h_sizes_rect, inhom_rect_file, ref_rel_diff, true, false);

        // Rectangle Inhomogeneous field with Material
        evaluate_jacobian_difference<inhom_field_rect_propagator_t>(
            track_count, euler_angles_I, euler_angles_F, rect_det_w_mat,
            detector_length, rect_bparam, inhom_bfield, overstep_tol,
            on_surface_tol, rk_tol_jac, rk_tol_dis, constraint_step_size,
            h_sizes_rect, inhom_rect_material_file, ref_rel_diff, true, false);

        // Rectangle Inhomogeneous field with Material (Covariance
        // transport)
        evaluate_covariance_transport<inhom_field_rect_propagator_t>(
            track_count, euler_angles_I, euler_angles_F, rect_det_w_mat,
            detector_length, rect_bparam, inhom_bfield, overstep_tol,
            on_surface_tol, rk_tol_cov, rk_tol_dis, constraint_step_size,
            rect_cov_transport_file, true);
      }
    }

    /**********************************
     * Wire telescope geometry
     **********************************/

    if (verbose_lvl >= 3 && !skip_wire) {
      std::clog << "Simulating wire telescope..." << std::endl;
    }

    // Get initial parameter
    const auto wire_bparam =
        get_initial_parameter<decltype(wire_det),
                              decltype(wire_det)::masks::id::e_drift_cell>(
            wire_det, track, B_z, helix_tol);

    if (!skip_wire) {
      scalar ref_rel_diff;

      if (rk_tolerance_iterate_mode) {
        std::array<unsigned int, 5u> num_iterations;
        std::array<bool, 25u> convergence;
        auto differentiated_jacobian =
            directly_differentiate<inhom_field_wire_propagator_t>(
                track_count, wire_bparam, wire_det_w_mat, detector_length,
                inhom_bfield, overstep_tol, on_surface_tol, rk_tol_dis,
                constraint_step_size, h_sizes_wire, num_iterations,
                convergence);

        for (std::size_t i = 0u; i < log10_tols.size(); i++) {
          // Wire Inhomogeneous field with Material
          evaluate_jacobian_difference<inhom_field_wire_propagator_t>(
              track_count, euler_angles_I, euler_angles_F, wire_det_w_mat,
              detector_length, wire_bparam, inhom_bfield, overstep_tol,
              on_surface_tol, std::pow(10.f, log10_tols[i]), rk_tol_dis,
              constraint_step_size, h_sizes_wire, wire_files[i], ref_rel_diff,
              true, do_inspect, true, differentiated_jacobian, num_iterations,
              convergence);

          dqopdqop_rel_diffs_wire[i].push_back(ref_rel_diff);
        }
      } else if (!rk_tolerance_iterate_mode) {
        // For helix
        evaluate_jacobian_difference_helix<
            decltype(wire_det), decltype(wire_det)::masks::id::e_drift_cell>(
            track_count, euler_angles_I, euler_angles_F, wire_det,
            detector_length, wire_bparam, B_z, h_sizes_wire, helix_wire_file,
            helix_tol);

        // Wire Const field
        evaluate_jacobian_difference<const_field_wire_propagator_t>(
            track_count, euler_angles_I, euler_angles_F, wire_det,
            detector_length, wire_bparam, const_bfield, overstep_tol,
            on_surface_tol, rk_tol_jac, rk_tol_dis, constraint_step_size,
            h_sizes_wire, const_wire_file, ref_rel_diff, true, false);

        // Wire Inhomogeneous field
        evaluate_jacobian_difference<inhom_field_wire_propagator_t>(
            track_count, euler_angles_I, euler_angles_F, wire_det,
            detector_length, wire_bparam, inhom_bfield, overstep_tol,
            on_surface_tol, rk_tol_jac, rk_tol_dis, constraint_step_size,
            h_sizes_wire, inhom_wire_file, ref_rel_diff, true, false);

        // Wire Inhomogeneous field with Material
        evaluate_jacobian_difference<inhom_field_wire_propagator_t>(
            track_count, euler_angles_I, euler_angles_F, wire_det_w_mat,
            detector_length, wire_bparam, inhom_bfield, overstep_tol,
            on_surface_tol, rk_tol_jac, rk_tol_dis, constraint_step_size,
            h_sizes_wire, inhom_wire_material_file, ref_rel_diff, true, false);

        // Wire Inhomogeneous field with Material (Covariance transport)
        evaluate_covariance_transport<inhom_field_wire_propagator_t>(
            track_count, euler_angles_I, euler_angles_F, wire_det_w_mat,
            detector_length, wire_bparam, inhom_bfield, overstep_tol,
            on_surface_tol, rk_tol_cov, rk_tol_dis, constraint_step_size,
            wire_cov_transport_file, true);
      }
    }
  }

  if (rk_tolerance_iterate_mode) {
    for (std::size_t i = 0u; i < log10_tols.size(); i++) {
      EXPECT_EQ(dqopdqop_rel_diffs_rect[i].size(), n_tracks);
      EXPECT_EQ(dqopdqop_rel_diffs_wire[i].size(), n_tracks);

      EXPECT_GE(statistics::mean(dqopdqop_rel_diffs_rect[i]), 1e-12f);
      EXPECT_LE(statistics::mean(dqopdqop_rel_diffs_rect[i]), 1e-2f);
      EXPECT_GE(statistics::mean(dqopdqop_rel_diffs_wire[i]), 1e-12f);
      EXPECT_LE(statistics::mean(dqopdqop_rel_diffs_wire[i]), 1e-2f);
    }
  }

  // Close files
  helix_rect_file.close();
  const_rect_file.close();
  inhom_rect_file.close();
  inhom_rect_material_file.close();

  helix_wire_file.close();
  const_wire_file.close();
  inhom_wire_file.close();
  inhom_wire_material_file.close();

  rect_cov_transport_file.close();
  wire_cov_transport_file.close();

  if (rk_tolerance_iterate_mode) {
    for (std::size_t i = 0u; i < log10_tols.size(); i++) {
      rect_files[i].close();
      wire_files[i].close();
    }
  }
}

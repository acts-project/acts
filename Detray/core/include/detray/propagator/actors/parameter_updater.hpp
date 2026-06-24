// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/geometry/tracking_surface.hpp"
#include "detray/propagator/detail/codegen/covariance_transport.hpp"
#include "detray/propagator/detail/codegen/full_jacobian.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/composite_actor.hpp"
#include "detray/propagator/detail/jacobian_engine.hpp"
#include "detray/propagator/detail/noise_estimation.hpp"
#include "detray/propagator/propagation_config.hpp"
#include "detray/utils/curvilinear_frame.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/type_registry.hpp"

namespace detray::actor {

template <concepts::algebra algebra_t>
struct parameter_transporter;

template <concepts::algebra algebra_t>
struct parameter_setter;

/// Configuration of noise estimation
struct noise_cfg {
  /// Percentage of total track path to assume as accumulated error
  float accumulated_error{0.001f};
  /// Number of standard deviations to assume to model the scattering noise
  int n_stddev{2};
  /// Estimate mask tolerance for navigation to for compensate scattering
  bool estimate_scattering_noise{true};
};

/// State of the parameter updater. Used by two distinct steps: Transporter and
/// setter
template <concepts::algebra algebra_t>
struct parameter_updater_state {
  // Allow only the corresponding actors to modify the parameter data
  friend struct parameter_transporter<algebra_t>;
  friend struct parameter_setter<algebra_t>;

  constexpr parameter_updater_state() = default;

  /// Start without explicit track parameters
  DETRAY_HOST_DEVICE
  explicit constexpr parameter_updater_state(
      const propagation::config& cfg,
      bound_matrix<algebra_t>* full_jac = nullptr)
      : m_full_jacobian(full_jac),
        m_cfg{cfg.navigation.accumulated_error,
              cfg.navigation.n_scattering_stddev,
              cfg.navigation.estimate_scattering_noise} {}

  /// Start from free track parameters
  DETRAY_HOST_DEVICE
  constexpr parameter_updater_state(
      const propagation::config& cfg,
      const free_track_parameters<algebra_t>& free_params,
      bound_matrix<algebra_t>* full_jac = nullptr)
      : m_full_jacobian(full_jac),
        m_cfg{cfg.navigation.accumulated_error,
              cfg.navigation.n_scattering_stddev,
              cfg.navigation.estimate_scattering_noise} {
    // Set bound track parameters
    curvilinear_frame<algebra_t> cf(free_params);
    m_bound_params.set_parameter_vector(cf.m_bound_vec);
  }

  /// Start from bound track parameters
  DETRAY_HOST_DEVICE
  constexpr parameter_updater_state(
      const propagation::config& cfg,
      const bound_track_parameters<algebra_t>& bound_params,
      bound_matrix<algebra_t>* full_jac = nullptr)
      : m_bound_params{bound_params},
        m_full_jacobian(full_jac),
        m_cfg{cfg.navigation.accumulated_error,
              cfg.navigation.n_scattering_stddev,
              cfg.navigation.estimate_scattering_noise} {}

  /// Initialize the state from free track parameters
  DETRAY_HOST_DEVICE
  constexpr void init(const free_track_parameters<algebra_t>& free_params) {
    // Set bound track parameters
    curvilinear_frame<algebra_t> cf(free_params);
    m_bound_params.set_parameter_vector(cf.m_bound_vec);

    // A dummy covariance - should not be used
    m_bound_params.set_covariance(matrix::identity<bound_matrix<algebra_t>>());

    // An invalid geometry identifier - should not be used
    m_bound_params.set_surface_link(geometry::identifier{});

    assert(!m_bound_params.is_invalid());
  }

  /// Initialize the state from bound track parameters
  DETRAY_HOST_DEVICE
  constexpr void init(const bound_track_parameters<algebra_t>& bound_params) {
    assert(!bound_params.is_invalid());
    m_bound_params = bound_params;
  }

  /// Always update the track parameters in the stepper state. Otherwise,
  /// observing actors decide if the parameters should be updated.
  DETRAY_HOST_DEVICE
  void always_update(const bool do_update = true) {
    m_always_update = do_update;
  }

  /// Notify the observing actors on the initial propagation surface (i.e.
  /// before the stepper advanced the state from the seed parameters)
  DETRAY_HOST_DEVICE
  void notify_on_initial(const bool do_notify = true) {
    m_notify_on_initial = do_notify;
  }

  /// @return access to the noise estimation configuration - const
  DETRAY_HOST_DEVICE
  const noise_cfg& noise_estimation_cfg() const { return m_cfg; }

  /// @return access to the noise estimation configuration
  DETRAY_HOST_DEVICE
  noise_cfg& noise_estimation_cfg() { return m_cfg; }

  /// @returns bound track parameters - const
  /// @note Only meaningful, if the navigation is at the corresponding surface
  DETRAY_HOST_DEVICE
  constexpr const bound_track_parameters<algebra_t>& bound_params() const {
    return m_bound_params;
  }

  /// @returns bound track parameters
  /// @note Only meaningful, if the navigation is at the corresponding surface
  DETRAY_HOST_DEVICE
  constexpr bound_track_parameters<algebra_t>& bound_params() {
    return m_bound_params;
  }

  /// @returns true if the full Jacobian matrix should be assembled.
  DETRAY_HOST_DEVICE
  constexpr bool has_full_jacobian() const {
    return m_full_jacobian != nullptr;
  }

  /// @returns the current full Jacbian.
  DETRAY_HOST_DEVICE
  constexpr const bound_matrix<algebra_t>& full_jacobian() const {
    assert(has_full_jacobian());
    return *m_full_jacobian;
  }

  /// Set new full Jacbian.
  DETRAY_HOST_DEVICE
  constexpr void set_full_jacobian(const bound_matrix<algebra_t>& jac) {
    assert(has_full_jacobian());
    *m_full_jacobian = jac;
  }

  /// Set new full Jacbian.
  DETRAY_HOST_DEVICE
  constexpr void set_full_jacobian(bound_matrix<algebra_t>* jac_ptr) {
    assert(jac_ptr);
    m_full_jacobian = jac_ptr;
  }

 private:
  /// Bound track parameters and covariance
  bound_track_parameters<algebra_t> m_bound_params{};
  /// Full jacobian for up to the current destination surface
  bound_matrix<algebra_t>* m_full_jacobian{nullptr};
  /// Configuration for the noise estimation
  noise_cfg m_cfg{};
  /// Always update the free track parameters, even if nothing changed
  bool m_always_update{false};
  /// Notify observing actors on initial surface (propagation init)
  bool m_notify_on_initial{true};
};

/// Result of the param. transporter: bound track parameters at dest. sf.
template <concepts::algebra algebra_t>
struct parameter_transporter_result : public actor::result {
  constexpr parameter_transporter_result() = default;

  /// Parameterized constructor using a pointer to the destination parameters
  DETRAY_HOST_DEVICE
  constexpr parameter_transporter_result(
      actor::status stat, bound_track_parameters<algebra_t>* params_ptr)
      : actor::result{stat}, m_destination_params_ptr{params_ptr} {
    assert(m_destination_params_ptr);
  }

  /// Parameterized constructor using destination parameters
  DETRAY_HOST_DEVICE
  constexpr parameter_transporter_result(
      actor::status stat, bound_track_parameters<algebra_t>& params)
      : actor::result{stat}, m_destination_params_ptr{&params} {
    assert(m_destination_params_ptr);
  }

  /// @returns access to the destination bound track parameters - const
  DETRAY_HOST_DEVICE
  constexpr const bound_track_parameters<algebra_t>& destination_params()
      const {
    assert(this->status == actor::status::e_unknown ||
           m_destination_params_ptr);
    return *m_destination_params_ptr;
  }

  /// @returns access to the destination bound track parameters
  DETRAY_HOST_DEVICE
  constexpr bound_track_parameters<algebra_t>& destination_params() {
    assert(this->status == actor::status::e_unknown ||
           m_destination_params_ptr);
    return *m_destination_params_ptr;
  }

 private:
  /// @returns a string stream that prints the transporter result details
  DETRAY_HOST
  friend std::ostream& operator<<(
      std::ostream& os, const parameter_transporter_result<algebra_t>& res) {
    os << static_cast<actor::result>(res) << std::endl;
    os << "destination params:\n" << res.destination_params() << std::endl;
    return os;
  }

  /// Bound track parameters of destination surface
  bound_track_parameters<algebra_t>* m_destination_params_ptr{nullptr};
};

/// Transport the free track parameters at the destination surface
/// and the covariance at the departure surface to bound track parameters
/// at the destination surface (current sensitive/material surface)
template <concepts::algebra algebra_t>
struct parameter_transporter : base_actor {
  /// @name Type definitions for the struct
  /// @{
  // Transformation matching this struct
  using transform3_type = dtransform3D<algebra_t>;
  // The current free track parameters in the stepper algorithm
  using free_track_parameters_type = free_track_parameters<algebra_t>;
  // The track parameters bound to the current sensitive/material surface
  using bound_track_parameters_type = bound_track_parameters<algebra_t>;
  // Bound matrix type (bound covariance)
  using bound_matrix_type = bound_matrix<algebra_t>;
  /// @}

  /// Use the parameter updater state
  using state = parameter_updater_state<algebra_t>;
  using result = parameter_transporter_result<algebra_t>;

  /// Filter the masks of a detector according to the local frame type
  struct select_frame {
    template <typename mask_t>
    using type = typename mask_t::local_frame;
  };

  /// Visitors to the surface local coordinate frame
  /// @{
  struct get_bound_to_free_dpos_dloc_visitor {
    template <typename frame_t>
    DETRAY_HOST_DEVICE constexpr dmatrix<algebra_t, 3, 2> operator()(
        const frame_t& /*frame*/, const transform3_type& trf3,
        const free_track_parameters_type& params) const {
      return detail::jacobian_engine<algebra_t>::
          template bound_to_free_jacobian_submatrix_dpos_dloc<frame_t>(
              trf3, params.pos(), params.dir());
    }
  };

  struct get_bound_to_free_dpos_dangle_visitor {
    template <typename frame_t>
    DETRAY_HOST_DEVICE constexpr dmatrix<algebra_t, 3, 2> operator()(
        const frame_t& /*frame*/, const transform3_type& trf3,
        const free_track_parameters_type& params,
        const dmatrix<algebra_t, 3, 2>& ddir_dangle) const {
      return detail::jacobian_engine<algebra_t>::
          template bound_to_free_jacobian_submatrix_dpos_dangle<frame_t>(
              trf3, params.pos(), params.dir(), ddir_dangle);
    }
  };

  struct get_free_to_bound_dloc_dpos_visitor {
    template <typename frame_t, typename stepper_state_t>
    DETRAY_HOST_DEVICE constexpr dmatrix<algebra_t, 2, 3> operator()(
        const frame_t& /*frame*/, const transform3_type& trf3,
        const stepper_state_t& stepping) const {
      return detail::jacobian_engine<algebra_t>::
          template free_to_bound_jacobian_submatrix_dloc_dpos<frame_t>(
              trf3, stepping().pos(), stepping().dir());
    }
  };
  /// @}

  /// Actor interface
  template <typename propagator_state_t>
  DETRAY_HOST_DEVICE result operator()(state& updater_state,
                                       propagator_state_t& propagation) const {
    const auto& navigation = propagation.navigation();

    // Do covariance transport only when the track is on surface
    if (!(navigation.is_on_sensitive() ||
          navigation.encountered_sf_material())) {
      return {};
    }
    assert(navigation.is_on_surface());

    // Current track position and transport jacobian
    auto& stepping = propagation.stepping();

    // Destination surface (current)
    const auto dest_sf = navigation.current_surface();

    // Bound track params of departure surface (not yet updated)
    auto& departure_params = updater_state.bound_params();

    // Result to be passed on to observing actors
    result res{actor::status::e_notify, departure_params};

    // Covariance is transported only when the departure surface is an
    // actual tracking surface (i.e. this disables the covariance
    // transport from curvilinear frames).
    if (!departure_params.surface_link().is_invalid()) {
      DETRAY_DEBUG_HOST(
          "Actor: Departure surface: " << departure_params.surface_link());

      // There is no need to transport the identical/initial parameters
      assert(!dest_sf.identifier().is_invalid());
      if (departure_params.surface_link() == dest_sf.identifier()) {
        DETRAY_VERBOSE_HOST_DEVICE(
            "Actor: On initial surface (%d), parameter transport not "
            "required",
            dest_sf.index());

        // Notify observers?
        if (!updater_state.m_notify_on_initial) {
          res.status = actor::status::e_unknown;
        }

        return res;
      }

      DETRAY_DEBUG_HOST("Actor: Destination surface: " << dest_sf.identifier());
      DETRAY_VERBOSE_HOST_DEVICE(
          "Actor: Transport track covariance to surface %d", dest_sf.index());

      // Transport the covariance
      const bound_matrix<algebra_t> propagation_step_jacobian =
          get_full_jacobian(propagation, departure_params);

      // Update the full Jacobian, if required
      if (math::fabs(stepping.path_length()) > 0.f) {
        if (updater_state.has_full_jacobian()) {
          updater_state.set_full_jacobian(propagation_step_jacobian *
                                          updater_state.full_jacobian());
        }

        // Reset transport Jacobian to identity matrix
        stepping.reset_transport_jacobian();
      }

      // Copy old covariance in order to update new covariance in place
      const bound_matrix_type old_cov = departure_params.covariance();
      bound_matrix_type& new_cov = departure_params.covariance();

      detray::detail::transport_covariance_to_bound_impl(
          old_cov, propagation_step_jacobian, new_cov);
    } else {
      DETRAY_VERBOSE_HOST_DEVICE(
          "Actor: Departure surface link invalid, setting covariance to "
          "identity");

      // Dummy covariance that allows to multiply the transport jacobian
      departure_params.set_covariance(matrix::identity<bound_matrix_type>());
    }

    DETRAY_VERBOSE_HOST_DEVICE(
        "Actor: Convert track parameter vector to local for surface %d",
        dest_sf.index());

    // Get the bound parameters at the destination surface
    res.destination_params().set_parameter_vector(
        dest_sf.free_to_bound_vector(propagation.context(), stepping()));

    // Set new surface link
    res.destination_params().set_surface_link(dest_sf.identifier());

    assert(!res.destination_params().is_invalid());

    DETRAY_DEBUG_HOST("Actor: Transported bound param.:\n"
                      << res.destination_params());

    return res;
  }

  /// @returns the full Jacobian between departure and destination surfaces
  template <typename propagator_state_t>
  DETRAY_HOST_DEVICE constexpr bound_matrix_type get_full_jacobian(
      propagator_state_t& propagation,
      const bound_track_parameters_type& departure_params) const {
    // Map the surface shapes of the detector down to the common frames
    using detector_t = typename propagator_state_t::detector_type;
    using frame_registry_t =
        types::mapped_registry<typename detector_t::masks, select_frame>;

    const auto& gctx = propagation.context();
    const auto& stepping = propagation.stepping();
    const auto& navigation = propagation.navigation();

    // Our goal here is to compute the full Jacobian, which is given as:
    //
    // J_full = J_F2B * (D + I) * J_transport * J_B2F
    //
    // The transport Jacobian is given, but the free-to-bound and
    // bound-to-free matrices as well as the derivative matrix $D$ need
    // to be computed still. In order to avoid full-rank matrix
    // multiplications, we use the known substructure of the
    // aforementioned matrices in order to perform fewer operations. We
    // describe in detail the mathematics performed throughout the
    // remainder of this function.

    // Departure surface
    const tracking_surface dep_sf{navigation.detector(),
                                  departure_params.surface_link()};

    // Destination Surface
    const tracking_surface dest_sf = navigation.current_surface();

    // Free track params of departure surface
    const free_track_parameters<algebra_t> dep_free_params =
        dep_sf.bound_to_free_vector(gctx, departure_params);

    // Free track params of destination surface
    const free_track_parameters<algebra_t>& dest_free_params = stepping();

    // First, we will compute the bound-to-free Jacobian, which is given
    // by three sub-Jacobians. In particular, the matrix has an 8x6
    // shape and looks like this:
    //
    //        l0  l1 phi  th q/p   t
    //  px [[  A,  A,  B,  B,  0,  0],
    //  py  [  A,  A,  B,  B,  0,  0],
    //  pz  [  A,  A,  B,  B,  0,  0],
    //   t  [  0,  0,  0,  0,  0,  1],
    //  dx  [  0,  0,  C,  C,  0,  0],
    //  dy  [  0,  0,  C,  C,  0,  0],
    //  dz  [  0,  0,  0,  C,  0,  0],
    // q/p  [  0,  0,  0,  0,  1,  0]]
    //
    // Note that we thus only have 17 out of 48 non-trivial matrix
    // elements.
    //
    // In this matrix, A represents the d(pos)/d(loc) submatrix, B is
    // the d(pos)/d(angle) submatrix, and C is the d(dir)/d(angle)
    // submatrix, all of which are 3x2 in size. Also, A and B depend on
    // the frame type while C is computed in the same way for all
    // frames. Finally, note that submatrix B is the zero matrix for
    // most frame types.
    const auto& dep_trf3 = dep_sf.transform(gctx);
    const dmatrix<algebra_t, 3, 2> b2f_dpos_dloc =
        types::visit<frame_registry_t, get_bound_to_free_dpos_dloc_visitor>(
            dep_sf.shape_id(), dep_trf3, dep_free_params);

    const dmatrix<algebra_t, 3, 2> b2f_ddir_dangle =
        detail::jacobian_engine<algebra_t>::
            bound_to_free_jacobian_submatrix_ddir_dangle(departure_params);

    const dmatrix<algebra_t, 3, 2> b2f_dpos_dangle =
        types::visit<frame_registry_t, get_bound_to_free_dpos_dangle_visitor>(
            dep_sf.shape_id(), dep_trf3, dep_free_params, b2f_ddir_dangle);

    // Next, we compute the derivative which is defined as the outer
    // product of the two 8x1 vectors representing the path to free
    // derivative and the free to path derivative. Thus, it is the
    // product of the following two vectors:
    //
    //  px [[ A],  [[ B],^T
    //  py  [ A],   [ B],
    //  pz  [ A],   [ B],
    //   t  [ 0]]   [ 0],
    //  dx  [ A],   [ C],
    //  dy  [ A],   [ C],
    //  dz  [ A],   [ C],
    // q/p  [ A],   [ 0]]
    //
    // Where A is frame independent and non-zero, B is frame-dependent
    // and non-zero, while C is frame-dependent and non-zero only for
    // line frames.
    const dpoint3D<algebra_t> dest_glob_pos{dest_free_params.pos()};
    const dvector3D<algebra_t> dest_glob_dir{dest_free_params.dir()};

    auto vol = navigation.current_volume();
    const auto vol_mat_ptr =
        vol.has_material() ? vol.material_parameters(dest_glob_pos) : nullptr;

    const auto path_to_free_derivative =
        detail::jacobian_engine<algebra_t>::path_to_free_derivative(
            dest_glob_dir, stepping.dtds(), stepping.dqopds(vol_mat_ptr));

    const auto free_to_path_derivative = dest_sf.free_to_path_derivative(
        gctx, dest_glob_pos, dest_glob_dir, stepping.dtds());

    // Now, we compute the free-to-bound Jacobian which is of size 6x8
    // and has the following structure:
    //
    //        px  py  pz   t  dx  dy  dz q/p
    //  l0 [[  A,  A,  A,  0,  0,  0,  0,  0],
    //  l1  [  A,  A,  A,  0,  0,  0,  0,  0],
    // phi  [  0,  0,  0,  0,  B,  B,  0,  0],
    //  th  [  0,  0,  0,  0,  B,  B,  B,  0],
    // q/p  [  0,  0,  0,  0,  0,  0,  0,  1],
    //   t  [  0,  0,  0,  1,  0,  0,  0,  0]]
    //
    // Thus, the number of non-trivial elements is only 11 out of the 48
    // matrix elements. Also, submatrix A depends on the frame type
    // while submatrix B is the same for all frame types.
    const dmatrix<algebra_t, 2, 3> f2b_dloc_dpos =
        types::visit<frame_registry_t, get_free_to_bound_dloc_dpos_visitor>(
            dest_sf.shape_id(), dest_sf.transform(gctx), stepping);

    const dmatrix<algebra_t, 2, 3> f2b_dangle_ddir = detail::jacobian_engine<
        algebra_t>::free_to_bound_jacobian_submatrix_dangle_ddir(dest_glob_dir);

    // Finally, we can use our Sympy-generated full Jacobian computation
    // and return its result.
    bound_matrix_type full_jacobian;

    if constexpr (std::decay_t<propagator_state_t>::stepper_uses_gradient) {
      detail::update_full_jacobian_with_gradient_impl(
          stepping.internal_transport_jacobian(), b2f_dpos_dloc,
          b2f_ddir_dangle, b2f_dpos_dangle, path_to_free_derivative,
          free_to_path_derivative, f2b_dloc_dpos, f2b_dangle_ddir,
          full_jacobian);
    } else {
      detail::update_full_jacobian_without_gradient_impl(
          stepping.internal_transport_jacobian(), b2f_dpos_dloc,
          b2f_ddir_dangle, b2f_dpos_dangle, path_to_free_derivative,
          free_to_path_derivative, f2b_dloc_dpos, f2b_dangle_ddir,
          full_jacobian);
    }

    return full_jacobian;
  }
};

/// Set the free track parameters in the stepper state, in case
/// observing actors (e.g. the material interactor) changed the bound parameters
template <concepts::algebra algebra_t>
struct parameter_setter : base_actor {
  /// Access the same parameter updater state as the parameter transporter
  using state = parameter_updater_state<algebra_t>;

  /// Actor interface: Observer to the parameter transporter and its observers
  template <typename propagator_state_t>
  DETRAY_HOST_DEVICE void operator()(
      state& updater_state, propagator_state_t& propagation,
      const parameter_transporter_result<algebra_t>& res) const {
    // Update the track parameters only if necessary
    if (res.status != actor::status::e_success &&
        !updater_state.m_always_update) {
      return;
    }

    const auto& navigation = propagation.navigation();
    auto& stepping = propagation.stepping();

    // One of the transporter observers might have exited the navigation
    assert(navigation.is_on_surface() || !navigation.is_alive());

    DETRAY_VERBOSE_HOST_DEVICE("Actor: Update the track parameters");

    // Updated bound track parameters
    const bound_track_parameters<algebra_t>& bound_params =
        res.destination_params();

    // Update free params after bound params were changed by actors
    const tracking_surface dest_sf{navigation.detector(),
                                   bound_params.surface_link()};
    stepping() =
        dest_sf.bound_to_free_vector(propagation.context(), bound_params);

    assert(!stepping().is_invalid());

    // Track pos/dir is not always known precisely: adjust navigation
    // tolerances according to bound covariance
    if (const noise_cfg& cfg = updater_state.noise_estimation_cfg();
        cfg.estimate_scattering_noise) {
      detail::estimate_external_mask_tolerance(
          bound_params, propagation,
          static_cast<dscalar<algebra_t>>(cfg.n_stddev), cfg.accumulated_error);
    }

    DETRAY_DEBUG_HOST("Actor: Updated bound param.:\n" << bound_params);
  }
};

/// Call actors that depend on the bound track parameters safely together
/// with the parameter transporter and parameter setter
template <typename algebra_t, typename... transporter_observers>
using parameter_updater =
    composite_actor<parameter_transporter<algebra_t>, transporter_observers...,
                    parameter_setter<algebra_t>>;

}  // namespace detray::actor

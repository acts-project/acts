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
#include "detray/material/concepts.hpp"
#include "detray/material/detail/material_accessor.hpp"
#include "detray/navigation/caching_navigator.hpp"
#include "detray/navigation/detail/print_state.hpp"
#include "detray/propagator/actors.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/type_list.hpp"

// Detray IO include(s)
#include "detray/io/utils/create_path.hpp"
#include "detray/io/utils/file_handle.hpp"

// Test include(s)
#include "detray/test/utils/data_record.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <filesystem>

namespace detray::material_validator {

/// @brief Record the material budget per thickness or pathlength
template <concepts::scalar scalar_t>
struct track_material {
  /// Phi and eta values of the track for which the material was recorded
  /// @{
  scalar_t phi{detail::invalid_value<scalar_t>()};
  scalar_t eta{detail::invalid_value<scalar_t>()};
  /// @}
  /// Accumulated radiation length per pathlength through the material
  scalar_t sX0{0.f};
  /// Accumulated radiation length per thickness
  scalar_t tX0{0.f};
  /// Accumulated interaction length per pathlength through the material
  scalar_t sL0{0.f};
  /// Accumulated interaction length per thickness
  scalar_t tL0{0.f};
};

/// @brief Functor to retrieve the material parameters for a given local
/// position
struct get_material_record {
  template <typename mat_group_t, typename index_t, concepts::point2D point2_t,
            concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE auto operator()(
      [[maybe_unused]] const mat_group_t &mat_group,
      [[maybe_unused]] const index_t &index,
      [[maybe_unused]] const point2_t &loc,
      [[maybe_unused]] const scalar_t cos_inc_angle) const {
    using material_t = typename mat_group_t::value_type;

    constexpr auto inv{detail::invalid_value<scalar_t>()};
    constexpr geometry::identifier inv_sf{};

    // Access homogeneous surface material or material maps
    if constexpr (concepts::surface_material<material_t>) {
      // Slab or rod
      const auto mat = detail::material_accessor::get(mat_group, index, loc);

      // Empty material can occur in material maps, skip it
      if (!mat) {
        // Set the pathlength and thickness to zero so that they
        // are not counted
        return material_record<scalar_t>{inv_sf, 0.f, 0.f, inv, inv};
      }

      const scalar_t seg{mat.path_segment(cos_inc_angle, loc[0])};
      const scalar_t t{mat.thickness()};
      const scalar_t mat_X0{mat.get_material().X0()};
      const scalar_t mat_L0{mat.get_material().L0()};

      return material_record<scalar_t>{inv_sf, seg, t, mat_X0, mat_L0};
    } else {
      return material_record<scalar_t>{inv_sf, inv, inv, inv, inv};
    }
  }
};

/// @brief Actor that collects all material encountered by a track during
///        navigation
///
/// The material is scaled with either the slab thickness or pathlength through
/// the material.
template <concepts::scalar scalar_t, template <typename...> class vector_t>
struct material_tracer : public detray::base_actor {
  using track_material_type = track_material<scalar_t>;
  using material_record_type = material_record<scalar_t>;

  struct state {
    friend struct material_tracer;

    /// Construct the vector containers with a given resource
    /// @param resource
    DETRAY_HOST
    explicit state(vecmem::memory_resource &resource)
        : m_mat_steps(&resource) {}

    /// Construct from externally provided vector for the @param steps
    DETRAY_HOST_DEVICE
    explicit state(vector_t<material_record_type> &&steps)
        : m_mat_steps(std::move(steps)) {}

    /// Access to the total recorded material along the track - const
    DETRAY_HOST_DEVICE
    const auto &get_track_material() const { return m_track_material; }

    /// Move the total recorded material out of the actor
    DETRAY_HOST
    auto &&release_track_material() && { return std::move(m_track_material); }

    /// Access to the recorded material steps along the track - const
    DETRAY_HOST_DEVICE
    const auto &get_material_steps() const { return m_mat_steps; }

    /// Move the recorded material steps pout of the actor
    DETRAY_HOST_DEVICE
    auto &&release_material_steps() && { return std::move(m_mat_steps); }

   private:
    /// Accumulated material data for the track
    track_material_type m_track_material{};

    /// Collect material parameters for every step
    vector_t<material_record_type> m_mat_steps{};
  };

  /// Run as observer to the parameter transporter (covariance transport)
  template <typename propagator_state_t, concepts::algebra algebra_t>
  DETRAY_HOST_DEVICE void operator()(
      state &tracer, const propagator_state_t &prop_state,
      const detray::actor::parameter_transporter_result<algebra_t> &res) const {
    record_track_dir(tracer, prop_state);

    // Only count material if navigator encountered it
    const auto &navigation = prop_state.navigation();
    if (!navigation.encountered_sf_material()) {
      return;
    }

    // For now use default context
    typename propagator_state_t::detector_type::geometry_context gctx{};

    const auto &track_param = res.destination_params();
    dpoint2D<algebra_t> loc_pos = track_param.bound_local();

    record_mat_step(tracer, gctx, navigation.current_surface(), loc_pos,
                    track_param.dir());
  }

  /// Run in a propagation chain without covariance transport
  template <typename propagator_state_t>
  DETRAY_HOST_DEVICE void operator()(
      state &tracer, const propagator_state_t &prop_state) const {
    using algebra_t = typename propagator_state_t::detector_type::algebra_type;

    record_track_dir(tracer, prop_state);

    // Only count material if navigator encountered it
    const auto &navigation = prop_state.navigation();
    if (!navigation.encountered_sf_material()) {
      return;
    }

    // For now use default context
    typename propagator_state_t::detector_type::geometry_context gctx{};

    // Current surface
    const auto sf = navigation.current_surface();

    const auto &track_param = prop_state.stepping()();
    dvector3D<algebra_t> glob_dir = track_param.dir();
    dpoint2D<algebra_t> loc_pos =
        sf.global_to_bound(gctx, track_param.pos(), glob_dir);

    record_mat_step(tracer, gctx, sf, loc_pos, glob_dir);
  }

 private:
  /// Record the track direction
  template <typename propagator_state_t>
  DETRAY_HOST_DEVICE inline auto record_track_dir(
      state &tracer, const propagator_state_t &prop_state) const {
    using algebra_t = typename propagator_state_t::detector_type::algebra_type;
    using vector3_t = dvector3D<algebra_t>;

    // Record the initial track direction
    vector3_t glob_dir = prop_state.stepping()().dir();
    if (detray::detail::is_invalid_value(tracer.m_track_material.eta) &&
        detray::detail::is_invalid_value(tracer.m_track_material.phi)) {
      tracer.m_track_material.eta = vector::eta(glob_dir);
      tracer.m_track_material.phi = vector::phi(glob_dir);
    }
  }

  /// Record the data for a material step
  template <typename detector_t>
  DETRAY_HOST_DEVICE inline auto record_mat_step(
      state &tracer, const typename detector_t::geometry_context &gctx,
      const tracking_surface<detector_t> sf,
      const dpoint2D<typename detector_t::algebra_type> &loc_pos,
      const dvector3D<typename detector_t::algebra_type> &glob_dir) const {
    // Fetch the material parameters and pathlength through the material
    const auto mat_params = sf.template visit_material<get_material_record>(
        loc_pos, cos_angle(gctx, sf, glob_dir, loc_pos));

    const scalar_t seg{mat_params.path};
    const scalar_t t{mat_params.thickness};
    const scalar_t mx0{mat_params.mat_X0};
    const scalar_t ml0{mat_params.mat_L0};

    // Fill the material record
    if (mx0 > 0.f) {
      tracer.m_track_material.sX0 += seg / mx0;
      tracer.m_track_material.tX0 += t / mx0;
    }
    if (ml0 > 0.f) {
      tracer.m_track_material.sL0 += seg / ml0;
      tracer.m_track_material.tL0 += t / ml0;
    }
    if (t > 0.f) {
      tracer.m_mat_steps.push_back({sf.identifier(), seg, t, mx0, ml0});
    }
  }
};

/// Run the propagation and record test data along the way
template <typename detector_t>
inline auto record_material(
    const typename detector_t::geometry_context /*gctx*/,
    vecmem::memory_resource *host_mr, const detector_t &det,
    const propagation::config &cfg,
    const free_track_parameters<typename detector_t::algebra_type> &track) {
  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;

  using stepper_t = line_stepper<algebra_t>;
  using navigator_t = caching_navigator<detector_t>;

  // Propagator with pathlimit aborter
  using material_tracer_t =
      material_validator::material_tracer<scalar_t, vecmem::vector>;
  using pathlimit_aborter_t = actor::pathlimit_aborter<scalar_t>;
  using parameter_updater_t =
      actor::parameter_updater<algebra_t,
                               actor::pointwise_material_interactor<algebra_t>,
                               material_tracer_t>;
  using actor_chain_t = actor_chain<pathlimit_aborter_t, parameter_updater_t>;
  using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

  // Propagator
  propagator_t prop{cfg};

  // Build actor and propagator states
  typename pathlimit_aborter_t::state pathlimit_aborter_state{
      cfg.stepping.path_limit};
  actor::parameter_updater_state<algebra_t> updater_state{cfg};
  typename actor::pointwise_material_interactor<algebra_t>::state
      interactor_state{};
  typename material_tracer_t::state mat_tracer_state{*host_mr};

  auto actor_states = detray::tie(pathlimit_aborter_state, updater_state,
                                  interactor_state, mat_tracer_state);

  typename propagator_t::state propagation{track, det, cfg.context};

  // Run the propagation
  bool success = prop.propagate(propagation, actor_states);

  return std::make_tuple(success,
                         std::move(mat_tracer_state).release_track_material(),
                         std::move(mat_tracer_state).release_material_steps());
}

/// Compare the result of two material tracers
///
/// @param reference the reference trace to compare against
/// @param ref_record the reference total material to compare against (s/X0)
/// @param mat_trace the recoded material trace
/// @param mat_record the recoded total material along the track (s/X0)
/// @param trk_i the track index/identifier
/// @param rel_tol the maximum allowed relative error for any parameters
/// @param verbose debug output
///
/// @returns if the traces contains matching material steps and if matching
/// traces contain the same total material
template <concepts::scalar scalar_t>
inline auto compare_traces(
    const dvector<material_record<scalar_t>> &reference,
    const track_material<scalar_t> &ref_record,
    const dvector<material_record<scalar_t>> &mat_trace,
    const track_material<scalar_t> &mat_record,
    std::size_t trk_i = detail::invalid_value<std::size_t>(),
    const double rel_tol = 0.01, const bool verbose = true) {
  // Material traces contain records of different surfaces/material
  bool is_bad_comp{reference.size() != mat_trace.size()};
  // Overall recorded material is different (in case small errors accumulate)
  bool is_diff_mat{false};

  std::stringstream debug_msg{};
  debug_msg << "Track No. " << trk_i << ":\n----------------" << std::endl;

  if (is_bad_comp) {
    debug_msg << "-> Different no. of surfaces: " << mat_trace.size()
              << " (ref.: " << reference.size() << ")\n"
              << std::endl;
  } else {
    for (std::size_t j = 0u; j < math::min(reference.size(), mat_trace.size());
         ++j) {
      if (reference[j].geo_id != mat_trace[j].geo_id) {
        is_bad_comp = true;
        debug_msg << "-> Surfaces don't match: " << mat_trace[j].geo_id
                  << " (ref.: " << reference[j].geo_id << ")" << std::endl;
        continue;
      }

      // TODO: Use approx_equal from algebra-plugins
      // Compare thickness of the surface material
      if ((reference[j].thickness - mat_trace[j].thickness) /
              reference[j].thickness >
          rel_tol) {
        is_bad_comp = true;
        debug_msg << "-> On surface: " << reference[j].geo_id << ":"
                  << std::endl;
        debug_msg << "-> thickness: " << mat_trace[j].thickness
                  << ", thickness ref.: " << reference[j].thickness
                  << std::endl;
      }

      // Compare radiation length of the surface material
      if ((reference[j].mat_X0 - mat_trace[j].mat_X0) / reference[j].mat_X0 >
          rel_tol) {
        is_bad_comp = true;
        debug_msg << "-> On surface: " << reference[j].geo_id << ":"
                  << std::endl;
        debug_msg << "-> X0: " << mat_trace[j].mat_X0
                  << ", X0 ref.: " << reference[j].mat_X0 << std::endl;
      }

      // Compare interaction length of the surface material
      if ((reference[j].mat_L0 - mat_trace[j].mat_L0) / reference[j].mat_L0 >
          rel_tol) {
        is_bad_comp = true;
        debug_msg << "-> On surface: " << reference[j].geo_id << ":"
                  << std::endl;
        debug_msg << "-> L0: " << mat_trace[j].mat_L0
                  << ", L0 ref.: " << reference[j].mat_L0 << std::endl;
      }

      // Compare path of the track through the surface material
      if ((reference[j].path - mat_trace[j].path) / reference[j].path >
          rel_tol) {
        is_bad_comp = true;
        debug_msg << "-> On surface: " << reference[j].geo_id << ":"
                  << std::endl;
        debug_msg << "-> Mat. path: " << mat_trace[j].path / unit<scalar_t>::mm
                  << " mm, mat. path ref.: "
                  << reference[j].path / unit<scalar_t>::mm << " mm"
                  << std::endl;
      }
    }
  }

  // Compare the total accumulated material along the track
  const double ref_mat_X0{static_cast<double>(ref_record.sX0)};
  const double mat_X0{static_cast<double>(mat_record.sX0)};
  const double rel_error{(ref_mat_X0 - mat_X0) / ref_mat_X0};
  const bool small_mat{ref_mat_X0 < rel_tol && mat_X0 < rel_tol};

  // If almost no material was collected, the relative error can be
  // large, but also has little consequence for tracking
  if (!(small_mat || (rel_error <= rel_tol))) {
    // Already know that material cannot add up: Only flag the additional
    // cases here
    if (!is_bad_comp) {
      is_diff_mat = true;
    }
    debug_msg << "\nTotal material discrepancy of " << 100. * rel_error << "%\n"
              << std::endl;
  }

  if (verbose && (is_bad_comp || is_diff_mat)) {
    DETRAY_INFO_HOST(debug_msg.str());
  }

  return std::make_tuple(is_bad_comp, is_diff_mat);
}

/// Write the accumulated material of a track from @param mat_records to a csv
/// file to the path @param mat_file_name
template <concepts::scalar scalar_t>
auto write_material(const std::string &mat_file_name,
                    const dvector<track_material<scalar_t>> &mat_records) {
  const auto file_path = std::filesystem::path{mat_file_name};
  assert(file_path.extension() == ".csv");

  // Make sure path to file exists
  io::create_path(file_path.parent_path());

  detray::io::file_handle outfile{
      mat_file_name, std::ios::out | std::ios::binary | std::ios::trunc};
  *outfile << "eta,phi,mat_sX0,mat_sL0,mat_tX0,mat_tL0" << std::endl;

  for (const auto &rec : mat_records) {
    *outfile << rec.eta << "," << rec.phi << "," << rec.sX0 << "," << rec.sL0
             << "," << rec.tX0 << "," << rec.tL0 << std::endl;
  }
}

}  // namespace detray::material_validator

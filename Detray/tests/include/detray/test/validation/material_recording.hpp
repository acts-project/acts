// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/core/concepts.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/material/concepts.hpp"
#include "detray/material/detail/material_accessor.hpp"
#include "detray/propagator/actors.hpp"
#include "detray/utils/logging.hpp"

// Detray IO include(s)
#include "detray/io/utils/create_path.hpp"
#include "detray/io/utils/file_handle.hpp"

// Test include(s)
#include "detray/test/utils/data_record.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <array>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>

namespace detray::validation {

/// @brief Functor to retrieve the material parameters for a given local
/// position
struct get_material_record {
  /// Retrieve the material record at a local position.
  /// @param mat_group Material group to inspect.
  /// @param index Material index.
  /// @param loc Local position.
  /// @param cos_inc_angle Cosine of the incidence angle.
  /// @return Material record for the position.
  template <typename mat_group_t, typename index_t, concepts::point2D point2_t,
            concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE auto operator()(
      [[maybe_unused]] const mat_group_t &mat_group,
      [[maybe_unused]] const index_t &index,
      [[maybe_unused]] const point2_t &loc,
      [[maybe_unused]] const scalar_t cos_inc_angle) const {
    using material_t = typename mat_group_t::value_type;

    constexpr auto inv{detail::invalid_value<double>()};
    constexpr geometry::identifier inv_sf{};

    // Access homogeneous surface material or material maps
    if constexpr (concepts::surface_material<material_t>) {
      // Slab or rod
      const auto mat = detail::material_accessor::get(mat_group, index, loc);

      // Empty material can occur in material maps, skip it
      if (!mat.is_valid()) {
        // Set the pathlength and thickness to zero so that they are not counted
        return material_record{inv_sf, 0.f, 0.f, inv, inv};
      }

      const double seg{mat.path_segment(cos_inc_angle, loc[0])};
      const double t{mat.thickness()};
      const double mat_X0{mat.get_material().X0()};
      const double mat_L0{mat.get_material().L0()};

      return material_record{inv_sf, seg, t, mat_X0, mat_L0};
    } else {
      return material_record{inv_sf, 0., 0., inv, inv};
    }
  }
};

/// @brief Actor that collects all material encountered by a track during
///        navigation
///
/// The material is scaled with either the slab thickness or pathlength through
/// the material.
template <template <typename...> class vector_t>
struct material_tracer : public detray::base_actor {
  /// Type storing accumulated track material.
  using track_material_type = track_material;
  /// Type storing one material step.
  using material_record_type = material_record;

  /// State used to collect material during propagation.
  struct state {
    friend struct material_tracer;

    /// Construct the vector containers with a given resource
    /// @param resource Memory resource for the material steps.
    DETRAY_HOST
    explicit state(vecmem::memory_resource &resource)
        : m_mat_steps(&resource) {}

    /// Construct from externally provided vector for the @p steps.
    /// @param steps Vector containing the material steps.
    DETRAY_HOST_DEVICE
    explicit state(vector_t<material_record_type> &&steps)
        : m_mat_steps(std::move(steps)) {}

    /// Access to the total recorded material along the track - const
    /// @return Recorded track material.
    DETRAY_HOST_DEVICE
    const auto &get_track_material() const { return m_track_material; }

    /// Move the total recorded material out of the actor
    /// @return The moved track material.
    DETRAY_HOST
    auto &&release_track_material() && { return std::move(m_track_material); }

    /// Access to the recorded material steps along the track - const
    /// @return Recorded material steps.
    DETRAY_HOST_DEVICE
    const auto &get_material_steps() const { return m_mat_steps; }

    /// Move the recorded material steps pout of the actor
    /// @return The moved material steps.
    DETRAY_HOST_DEVICE
    auto &&release_material_steps() && { return std::move(m_mat_steps); }

   private:
    /// Accumulated material data for the track
    track_material_type m_track_material{};

    /// Collect material parameters for every step
    vector_t<material_record_type> m_mat_steps{};
  };

  /// Run as observer to the parameter transporter (covariance transport)
  /// @param tracer State collecting material.
  /// @param prop_state Current propagator state.
  /// @param res Parameter transport result.
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
  /// @param tracer State collecting material.
  /// @param prop_state Current propagator state.
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
  template <concepts::detector detector_t>
  DETRAY_HOST_DEVICE inline auto record_mat_step(
      state &tracer, const typename detector_t::geometry_context &gctx,
      const tracking_surface<detector_t> sf,
      const dpoint2D<typename detector_t::algebra_type> &loc_pos,
      const dvector3D<typename detector_t::algebra_type> &glob_dir) const {
    // Fetch the material parameters and pathlength through the material
    const auto mat_params = sf.template visit_material<get_material_record>(
        loc_pos, cos_angle(gctx, sf, glob_dir, loc_pos));

    const double seg{mat_params.path};
    const double t{mat_params.thickness};
    const double mx0{mat_params.mat_X0};
    const double ml0{mat_params.mat_L0};

    // Fill the material record
    if (mx0 > 0.) {
      tracer.m_track_material.sX0 += seg / mx0;
      tracer.m_track_material.tX0 += t / mx0;
    }
    if (ml0 > 0.) {
      tracer.m_track_material.sL0 += seg / ml0;
      tracer.m_track_material.tL0 += t / ml0;
    }
    if (t > 0.) {
      tracer.m_mat_steps.push_back({sf.identifier(), seg, t, mx0, ml0});
    }
  }
};

/// Name of a material data file: "<detector>_<file stem>[_<suffix>].csv"
///
/// @param det_name the name of the detector
/// @param material_file the configured output file (only the parent path and
///                      the stem are used)
/// @param suffix optional suffix, e.g. the name of the compute backend
///
/// @returns the full path to the material data file
inline std::filesystem::path material_file_path(
    const std::string &det_name, const std::string &material_file,
    std::string_view suffix = {}) {
  const std::filesystem::path mat_path{material_file};

  std::string file_name{det_name + "_" + mat_path.stem().string()};
  if (!suffix.empty()) {
    file_name += "_" + std::string(suffix);
  }

  return mat_path.parent_path() / (file_name + ".csv");
}

/// Write the accumulated material of a track from @param mat_records to a csv
/// file to the path @param mat_file_name
inline auto write_track_material(const std::string &mat_file_name,
                                 const dvector<track_material> &mat_records) {
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

/// Read the accumulated material per track from a csv file that was written
/// by @c write_material
///
/// @param mat_file_name the path to the csv file
///
/// @returns the accumulated material of every track in the file
inline auto read_track_material(const std::string &mat_file_name) {
  detray::io::file_handle infile{mat_file_name, std::ios::in};

  dvector<track_material> mat_records{};

  // Skip the header
  std::string line{};
  std::getline(*infile, line);

  while (std::getline(*infile, line)) {
    if (line.empty()) {
      continue;
    }

    std::istringstream line_stream{line};
    std::array<double, 6u> values{};
    for (std::size_t i = 0u; i < values.size(); ++i) {
      std::string value{};
      if (!std::getline(line_stream, value, ',')) {
        throw std::invalid_argument("Invalid line in material file " +
                                    mat_file_name + ": " + line);
      }
      values[i] = std::stod(value);
    }

    track_material rec{};
    rec.eta = values[0];
    rec.phi = values[1];
    rec.sX0 = values[2];
    rec.sL0 = values[3];
    rec.tX0 = values[4];
    rec.tL0 = values[5];

    mat_records.push_back(rec);
  }

  return mat_records;
}

}  // namespace detray::validation

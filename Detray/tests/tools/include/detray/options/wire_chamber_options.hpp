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

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/geometry/shapes/line.hpp"

// Detray test include(s)
#include "detray/options/options_handling.hpp"
#include "detray/test/common/build_wire_chamber.hpp"

// Boost
#include "detray/options/boost_program_options.hpp"

// System include(s)
#include <string>

namespace detray::options {

namespace detail {

/// Add options that are independent of the wire surface shape
template <concepts::scalar scalar_t, typename wire_shape_t>
void add_wire_chamber_options(
    boost::program_options::options_description &desc,
    const wire_chamber_config<scalar_t, wire_shape_t> &cfg) {
  desc.add_options()(
      "layers",
      boost::program_options::value<unsigned int>()->default_value(
          cfg.n_layers()),
      "number of layers")("half_z",
                          boost::program_options::value<float>()->default_value(
                              static_cast<float>(cfg.half_z())),
                          "half length z of the chamber [mm]")(
      "cell_size",
      boost::program_options::value<float>()->default_value(
          static_cast<float>(cfg.cell_size())),
      "half length of the wire cells/straw tubes [mm]")(
      "stereo_angle",
      boost::program_options::value<float>()->default_value(
          static_cast<float>(cfg.stereo_angle())),
      "abs. stereo angle [rad]")(
      "mat_radius",
      boost::program_options::value<float>()->default_value(
          static_cast<float>(cfg.mat_radius())),
      "radius of material rods [mm]");
}

/// Configure options that are independent of the wire surface shape
template <concepts::scalar scalar_t, typename wire_shape_t>
void configure_wire_chamber_options(
    const boost::program_options::variables_map &vm,
    wire_chamber_config<scalar_t, wire_shape_t> &cfg) {
  cfg.n_layers(vm["layers"].as<unsigned int>());
  cfg.half_z(vm["half_z"].as<float>());
  cfg.cell_size(vm["cell_size"].as<float>());
  cfg.stereo_angle(vm["stereo_angle"].as<float>());
  cfg.mat_radius(vm["mat_radius"].as<float>());
}

}  // namespace detail

/// Add options for the wire chamber with wire cells
/// @{
template <>
void add_options<wire_chamber_config<float, detray::line_square>>(
    boost::program_options::options_description &desc,
    const wire_chamber_config<float, detray::line_square> &cfg) {
  detail::add_wire_chamber_options(desc, cfg);
}
template <>
void add_options<wire_chamber_config<double, detray::line_square>>(
    boost::program_options::options_description &desc,
    const wire_chamber_config<double, detray::line_square> &cfg) {
  detail::add_wire_chamber_options(desc, cfg);
}
/// @}

/// Add options for the wire chamber with straw tubes
/// @{
template <>
void add_options<wire_chamber_config<float, detray::line_circular>>(
    boost::program_options::options_description &desc,
    const wire_chamber_config<float, detray::line_circular> &cfg) {
  detail::add_wire_chamber_options(desc, cfg);
}
template <>
void add_options<wire_chamber_config<double, detray::line_circular>>(
    boost::program_options::options_description &desc,
    const wire_chamber_config<double, detray::line_circular> &cfg) {
  detail::add_wire_chamber_options(desc, cfg);
}
/// @}

/// Configure the detray wire chamber with wire cells
/// @{
template <>
void configure_options<wire_chamber_config<float, detray::line_square>>(
    const boost::program_options::variables_map &vm,
    wire_chamber_config<float, detray::line_square> &cfg) {
  detail::configure_wire_chamber_options(vm, cfg);
}
template <>
void configure_options<wire_chamber_config<double, detray::line_square>>(
    const boost::program_options::variables_map &vm,
    wire_chamber_config<double, detray::line_square> &cfg) {
  detail::configure_wire_chamber_options(vm, cfg);
}
/// @}

/// Configure the detray wire chamber with straw tubes
/// @{
template <>
void configure_options<wire_chamber_config<float, detray::line_circular>>(
    const boost::program_options::variables_map &vm,
    wire_chamber_config<float, detray::line_circular> &cfg) {
  detail::configure_wire_chamber_options(vm, cfg);
}
template <>
void configure_options<wire_chamber_config<double, detray::line_circular>>(
    const boost::program_options::variables_map &vm,
    wire_chamber_config<double, detray::line_circular> &cfg) {
  detail::configure_wire_chamber_options(vm, cfg);
}
/// @}

}  // namespace detray::options

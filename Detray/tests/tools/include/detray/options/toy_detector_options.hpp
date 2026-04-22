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

// Detray test include(s)
#include "detray/options/options_handling.hpp"
#include "detray/test/common/build_toy_detector.hpp"

// Boost
#include "detray/options/boost_program_options.hpp"

// System include(s)
#include <stdexcept>
#include <string>

namespace detray::options {

namespace detail {

/// Add options for the detray toy detector
template <concepts::scalar scalar_t>
void add_toy_det_options(boost::program_options::options_description &desc,
                         const toy_det_config<scalar_t> &cfg) {
  desc.add_options()(
      "barrel_layers",
      boost::program_options::value<unsigned int>()->default_value(
          cfg.n_brl_layers()),
      "number of barrel layers [0-4]")(
      "endcap_layers",
      boost::program_options::value<unsigned int>()->default_value(
          cfg.n_edc_layers()),
      "number of endcap layers on either side [0-7]")(
      "homogeneous_material",
      "Generate homogeneous material description (default)")(
      "material_maps", "Generate material maps");
}

/// Configure the detray toy detector
template <concepts::scalar scalar_t>
void configure_toy_det_options(const boost::program_options::variables_map &vm,
                               toy_det_config<scalar_t> &cfg) {
  cfg.n_brl_layers(vm["barrel_layers"].as<unsigned int>());
  cfg.n_edc_layers(vm["endcap_layers"].as<unsigned int>());

  if (vm.count("homogeneous_material") && vm.count("material_maps")) {
    throw std::invalid_argument("Please specify only one material description");
  }
  if (vm.count("homogeneous_material")) {
    cfg.use_material_maps(false);
  }
  if (vm.count("material_maps")) {
    cfg.use_material_maps(true);
  }
}

}  // namespace detail

/// Add options for the toy detector
/// @{
template <>
void add_options<toy_det_config<float>>(
    boost::program_options::options_description &desc,
    const toy_det_config<float> &cfg) {
  detail::add_toy_det_options(desc, cfg);
}

template <>
void add_options<toy_det_config<double>>(
    boost::program_options::options_description &desc,
    const toy_det_config<double> &cfg) {
  detail::add_toy_det_options(desc, cfg);
}
/// @}

/// Configure the detray toy detector
/// @{
template <>
void configure_options<toy_det_config<float>>(
    const boost::program_options::variables_map &vm,
    toy_det_config<float> &cfg) {
  detail::configure_toy_det_options(vm, cfg);
}

template <>
void configure_options<toy_det_config<double>>(
    const boost::program_options::variables_map &vm,
    toy_det_config<double> &cfg) {
  detail::configure_toy_det_options(vm, cfg);
}
/// @}

}  // namespace detray::options

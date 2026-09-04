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
#include "detray/utils/logging.hpp"

// Detray test include(s)
#include "detray/options/options_handling.hpp"
#include "detray/test/common/event_generator/random_track_generator_config.hpp"
#include "detray/test/common/event_generator/uniform_track_generator_config.hpp"

// Boost
#include "detray/options/boost_program_options.hpp"

// System include(s)
#include <stdexcept>
#include <vector>

namespace detray::options {

namespace detail {

/// Add options for detray event generation
template <concepts::scalar scalar_t>
void add_uniform_track_gen_options(
    boost::program_options::options_description &desc,
    const uniform_track_generator_config<scalar_t> &cfg) {
  desc.add_options()(
      "phi_steps",
      boost::program_options::value<std::size_t>()->default_value(
          cfg.phi_steps()),
      "No. phi steps for particle gun")(
      "eta_steps",
      boost::program_options::value<std::size_t>()->default_value(
          cfg.eta_steps()),
      "No. eta steps for particle gun")(
      "random_seed",
      boost::program_options::value<std::size_t>()->default_value(cfg.seed()),
      "Seed for the random number generator")(
      "eta_range",
      boost::program_options::value<std::vector<scalar_t>>()->multitoken(),
      "Min, Max range of eta values for particle gun")(
      "randomize_charge", "Randomly flip charge sign per track")(
      "origin",
      boost::program_options::value<std::vector<scalar_t>>()->multitoken(),
      "Coordinates for particle gun origin position [mm]")(
      "p_range",
      boost::program_options::value<std::vector<scalar_t>>()->multitoken(),
      "Total momentum [range] of the test particles [GeV]")(
      "pT_range",
      boost::program_options::value<std::vector<scalar_t>>()->multitoken(),
      "Transverse momentum [range] of the test particles [GeV]");
}

/// Add options for detray event generation
template <concepts::scalar scalar_t>
void configure_uniform_track_gen_options(
    const boost::program_options::variables_map &vm,
    uniform_track_generator_config<scalar_t> &cfg) {
  cfg.phi_steps(vm["phi_steps"].as<std::size_t>());
  cfg.eta_steps(vm["eta_steps"].as<std::size_t>());
  cfg.seed(vm["random_seed"].as<std::size_t>());
  cfg.randomize_charge(vm.count("randomize_charge"));

  if (vm.count("eta_range") != 0u) {
    const auto eta_range = vm["eta_range"].as<std::vector<scalar_t>>();
    if (eta_range.size() == 2u) {
      cfg.eta_range(eta_range[0], eta_range[1]);
    } else {
      throw std::invalid_argument("Eta range needs two arguments");
    }
  }
  if (vm.count("origin") != 0u) {
    const auto origin = vm["origin"].as<std::vector<scalar_t>>();
    if (origin.size() == 3u) {
      cfg.origin(origin[0] * unit<scalar_t>::mm, origin[1] * unit<scalar_t>::mm,
                 origin[2] * unit<scalar_t>::mm);
    } else {
      throw std::invalid_argument("Particle gun origin needs three arguments");
    }
  }
  if (vm.count("pT_range") != 0u && vm.count("p_range") != 0u) {
    throw std::invalid_argument(
        "Transverse and total momentum cannot be specified at the same "
        "time");
  }
  if (vm.count("pT_range") != 0u) {
    const auto pt_range = vm["pT_range"].as<std::vector<scalar_t>>();

    // Default
    if (pt_range.empty()) {
      cfg.p_T(static_cast<scalar_t>(cfg.m_p_mag));
    } else if (pt_range.size() <= 2u) {
      if (pt_range.size() == 2u) {
        DETRAY_ERROR_HOST(
            "Momentum range not possible with uniform "
            "track generator: Using first value.");
      }
      cfg.p_T(pt_range[0] * unit<scalar_t>::GeV);
    } else {
      throw std::invalid_argument(
          "Wrong number of arguments for pT range: Need one argument or "
          "range");
    }
  } else {
    auto p_range = std::vector<scalar_t>{};
    if (vm.count("pT_range") != 0u) {
      p_range = vm["p_range"].as<std::vector<scalar_t>>();
    }

    // Default
    if (p_range.empty()) {
      cfg.p_tot(static_cast<scalar_t>(cfg.m_p_mag));
    } else if (p_range.size() <= 2u) {
      if (p_range.size() == 2u) {
        DETRAY_ERROR_HOST(
            "Momentum range not possible with uniform "
            "track generator: Using first value.");
      }
      cfg.p_tot(p_range[0] * unit<scalar_t>::GeV);
    } else {
      std::string err_str{
          "Wrong number of arguments for p_tot range: Need one argument "
          "or range"};
      DETRAY_FATAL_HOST(err_str);
      throw std::invalid_argument(err_str);
    }
  }
}

/// Add options for detray event generation
template <concepts::scalar scalar_t>
void add_rnd_track_gen_options(
    boost::program_options::options_description &desc,
    const random_track_generator_config<scalar_t> &cfg) {
  desc.add_options()(
      "n_tracks",
      boost::program_options::value<std::size_t>()->default_value(
          cfg.n_tracks()),
      "No. of tracks for particle gun")(
      "random_seed",
      boost::program_options::value<std::size_t>()->default_value(cfg.seed()),
      "Seed for the random number generator")(
      "theta_range",
      boost::program_options::value<std::vector<scalar_t>>()->multitoken(),
      "Min, Max range of theta values for particle gun. Interval in [0, pi)")(
      "eta_range",
      boost::program_options::value<std::vector<scalar_t>>()->multitoken(),
      "Min, Max range of eta values for particle gun")(
      "randomize_charge", "Randomly flip charge sign per track")(
      "origin",
      boost::program_options::value<std::vector<scalar_t>>()->multitoken(),
      "Coordinates for particle gun origin position [mm]")(
      "p_range",
      boost::program_options::value<std::vector<scalar_t>>()->multitoken(),
      "Total momentum [range] of the test particles [GeV]")(
      "pT_range",
      boost::program_options::value<std::vector<scalar_t>>()->multitoken(),
      "Transverse momentum [range] of the test particles [GeV]");
}

/// Add options for detray event generation
template <concepts::scalar scalar_t>
void configure_rnd_track_gen_options(
    const boost::program_options::variables_map &vm,
    random_track_generator_config<scalar_t> &cfg) {
  cfg.n_tracks(vm["n_tracks"].as<std::size_t>());
  cfg.seed(vm["random_seed"].as<std::size_t>());
  cfg.randomize_charge(vm.count("randomize_charge"));

  if (vm.count("eta_range") != 0u && vm.count("theta_range") != 0u) {
    throw std::invalid_argument(
        "Eta range and theta range cannot be specified at the same time");
  } else if (vm.count("eta_range") != 0u) {
    const auto eta_range = vm["eta_range"].as<std::vector<scalar_t>>();
    if (eta_range.size() == 2u) {
      scalar_t min_theta{2.f * std::atan(std::exp(-eta_range[0]))};
      scalar_t max_theta{2.f * std::atan(std::exp(-eta_range[1]))};

      // Wrap around
      if (min_theta > max_theta) {
        scalar_t tmp{min_theta};
        min_theta = max_theta;
        max_theta = tmp;
      }

      cfg.theta_range(min_theta, max_theta);
    } else {
      throw std::invalid_argument("Eta range needs two arguments");
    }
  } else if (vm.count("theta_range") != 0u) {
    const auto theta_range = vm["theta_range"].as<std::vector<scalar_t>>();
    if (theta_range.size() == 2u) {
      cfg.theta_range(theta_range[0], theta_range[1]);
    } else {
      throw std::invalid_argument("Theta range needs two arguments");
    }
  }
  if (vm.count("origin") != 0u) {
    const auto origin = vm["origin"].as<std::vector<scalar_t>>();
    if (origin.size() == 3u) {
      cfg.origin(origin[0] * unit<scalar_t>::mm, origin[1] * unit<scalar_t>::mm,
                 origin[2] * unit<scalar_t>::mm);
    } else {
      throw std::invalid_argument(
          "Particle gun origin needs three coordinates");
    }
  }
  if (vm.count("pT_range") != 0u && vm.count("p_range") != 0u) {
    throw std::invalid_argument(
        "Transverse and total momentum cannot be specified at the same "
        "time");
  }
  if (vm.count("pT_range") != 0u) {
    const auto pt_range = vm["pT_range"].as<std::vector<scalar_t>>();

    // Default
    if (pt_range.empty()) {
      cfg.p_T(cfg.mom_range()[0]);
    } else if (pt_range.size() == 1u) {
      cfg.p_T(pt_range[0] * unit<scalar_t>::GeV);
    } else if (pt_range.size() == 2u) {
      cfg.pT_range(pt_range[0] * unit<scalar_t>::GeV,
                   pt_range[1] * unit<scalar_t>::GeV);
    } else {
      throw std::invalid_argument(
          "Wrong number of arguments for pT range: Need one argument or "
          "range");
    }
  } else {
    auto p_range = std::vector<scalar_t>();
    if (vm.count("pT_range") != 0u) {
      p_range = vm["p_range"].as<std::vector<scalar_t>>();
    }

    // Default
    if (p_range.empty()) {
      cfg.p_tot(cfg.mom_range()[0]);
    } else if (p_range.size() == 1u) {
      cfg.p_tot(p_range[0] * unit<scalar_t>::GeV);
    } else if (p_range.size() == 2u) {
      cfg.mom_range(p_range[0] * unit<scalar_t>::GeV,
                    p_range[1] * unit<scalar_t>::GeV);
    } else {
      throw std::invalid_argument(
          "Wrong number of arguments for p_tot range: Need one argument "
          "or range");
    }
  }
}

}  // namespace detail

/// Add options for the uniform track generator
/// @{
template <>
void add_options<uniform_track_generator_config<float>>(
    boost::program_options::options_description &desc,
    const uniform_track_generator_config<float> &cfg) {
  detail::add_uniform_track_gen_options(desc, cfg);
}

template <>
void add_options<uniform_track_generator_config<double>>(
    boost::program_options::options_description &desc,
    const uniform_track_generator_config<double> &cfg) {
  detail::add_uniform_track_gen_options(desc, cfg);
}
/// @}

/// Configure the detray uniform track generator
/// @{
template <>
void configure_options<uniform_track_generator_config<float>>(
    const boost::program_options::variables_map &vm,
    uniform_track_generator_config<float> &cfg) {
  detail::configure_uniform_track_gen_options(vm, cfg);
}

template <>
void configure_options<uniform_track_generator_config<double>>(
    const boost::program_options::variables_map &vm,
    uniform_track_generator_config<double> &cfg) {
  detail::configure_uniform_track_gen_options(vm, cfg);
}
/// @}

/// Add options for the random track generator
/// @{
template <>
void add_options<random_track_generator_config<float>>(
    boost::program_options::options_description &desc,
    const random_track_generator_config<float> &cfg) {
  detail::add_rnd_track_gen_options(desc, cfg);
}

template <>
void add_options<random_track_generator_config<double>>(
    boost::program_options::options_description &desc,
    const random_track_generator_config<double> &cfg) {
  detail::add_rnd_track_gen_options(desc, cfg);
}
/// @}

/// Configure the detray random track generator
/// @{
template <>
void configure_options<random_track_generator_config<float>>(
    const boost::program_options::variables_map &vm,
    random_track_generator_config<float> &cfg) {
  detail::configure_rnd_track_gen_options(vm, cfg);
}

template <>
void configure_options<random_track_generator_config<double>>(
    const boost::program_options::variables_map &vm,
    random_track_generator_config<double> &cfg) {
  detail::configure_rnd_track_gen_options(vm, cfg);
}
/// @}

}  // namespace detray::options

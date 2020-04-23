// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/program_options.hpp>
#include <boost/progress.hpp>
#include <random>
#include <string>

#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "ACTFW/Utilities/Options.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Units.hpp"

/// The main executable
///
/// Creates an InterpolatedBFieldMap from a txt or csv file
/// It then tests random versus stepwise access with the
/// direct getField access and the cell.getField access
/// with cell caching

namespace po = boost::program_options;

using UniformDist = std::uniform_real_distribution<double>;
using RandomEngine = std::mt19937;

template <typename field_t, typename field_context_t>
void accessStepWise(field_t& bField, field_context_t& bFieldContext,
                    size_t events, size_t theta_steps, double theta_0,
                    double theta_step, size_t phi_steps, double phi_0,
                    double phi_step, size_t access_steps, double access_step) {
  std::cout << "[>>>] Start: step-wise access pattern ... " << std::endl;
  size_t mismatched = 0;
  // initialize the field cache
  typename field_t::Cache bCache(bFieldContext);
  // boost display
  size_t totalSteps = events * theta_steps * phi_steps * access_steps;
  boost::progress_display show_progress(totalSteps);
  // the event loop
  // loop over the events - @todo move to parallel for
  for (size_t ievt = 0; ievt < events; ++ievt) {
    for (size_t itheta = 0; itheta < theta_steps; ++itheta) {
      double theta = theta_0 + itheta * theta_step;
      for (size_t iphi = 0; iphi < phi_steps; ++iphi) {
        double phi = phi_0 + iphi * phi_step;
        // make a direction
        Acts::Vector3D dir(cos(phi) * sin(theta), sin(phi) * sin(theta),
                           cos(theta));
        // check for the current step
        double currentStep = 0.;
        // now step through the magnetic field
        for (size_t istep = 0; istep < access_steps; ++istep) {
          Acts::Vector3D position = currentStep * dir;
          // access the field directly
          auto field_direct = bField.getField(position);
          // access the field with the cell
          auto field_from_cache = bField.getField(position, bCache);
          // check
          if (!field_direct.isApprox(field_from_cache)) {
            ++mismatched;
          }
          // increase the step
          currentStep += access_step;
          // show the progress bar
          ++show_progress;
        }
      }
    }
    std::cout << "[<<<] End result : " << mismatched << "/" << totalSteps
              << " mismatches" << std::endl;
  }
}

template <typename field_t, typename field_context_t>
void accessRandom(field_t& bField, field_context_t& bFieldContext,
                  size_t totalSteps, double radius) {
  std::cout << "[>>>] Start: random access pattern ... " << std::endl;
  size_t mismatched = 0;
  RandomEngine rng;
  UniformDist xDist(-radius, radius);
  UniformDist yDist(-radius, radius);
  UniformDist zDist(-radius, radius);

  // initialize the field cache
  typename field_t::Cache bCache(bFieldContext);
  boost::progress_display show_progress(totalSteps);

  // the event loop
  // loop over the events - @todo move to parallel for
  for (size_t istep = 0; istep < totalSteps; ++istep) {
    Acts::Vector3D position(xDist(rng), yDist(rng), zDist(rng));
    // access the field directly
    auto field_direct = bField.getField(position);
    // access the field with the cell
    auto field_from_cache = bField.getField(position, bCache);
    // check
    if (!field_direct.isApprox(field_from_cache)) {
      ++mismatched;
    }
    // show the progress bar
    ++show_progress;
  }
  std::cout << "[<<<] End result : " << mismatched << "/" << totalSteps
            << " mismatches" << std::endl;
}

/// @brief main executable
///
/// @param argc The argument count
/// @param argv The argument list
int main(int argc, char* argv[]) {
  // Declare the supported program options.
  auto desc = FW::Options::makeDefaultOptions();
  FW::Options::addSequencerOptions(desc);
  FW::Options::addBFieldOptions(desc);
  desc.add_options()("bf-phi-range",
                     po::value<read_range>()->default_value({-M_PI, M_PI}),
                     "range in which the phi parameter is generated.")(
      "bf-theta-range", po::value<read_range>()->default_value({0., M_PI}),
      "range in which the eta parameter is generated.")(
      "bf-phisteps", po::value<size_t>()->default_value(1000),
      "number of steps for the phi parameter.")(
      "bf-thetasteps", po::value<size_t>()->default_value(100),
      "number of steps for the eta parameter.")(
      "bf-accesssteps", po::value<size_t>()->default_value(100),
      "number of steps for magnetic field access.")(
      "bf-tracklength", po::value<double>()->default_value(100.),
      "track length in [mm] magnetic field access.");
  auto vm = FW::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  // A test magnetic field context
  Acts::MagneticFieldContext magFieldContext = Acts::MagneticFieldContext();

  // TODO
  // Why does this need number-of-events? If it really does emulate
  // per-event access patterns this should be switched to a proper
  // Sequencer-based tool. Otherwise it should be removed.
  auto nEvents = FW::Options::readSequencerConfig(vm).events;
  auto bFieldVar = FW::Options::readBField(vm);

  // Get the phi and eta range
  auto phir = vm["bf-phi-range"].as<read_range>();
  auto thetar = vm["bf-theta-range"].as<read_range>();
  // Get the granularity
  size_t phi_steps = vm["bf-phisteps"].as<size_t>();
  size_t theta_steps = vm["bf-thetasteps"].as<size_t>();
  // The defaults
  size_t access_steps = vm["bf-accesssteps"].as<size_t>();
  double track_length = vm["bf-tracklength"].as<double>() * Acts::units::_mm;
  // sort the ranges - and prepare the access grid
  std::sort(phir.begin(), phir.end());
  std::sort(thetar.begin(), thetar.end());
  double phi_span = std::abs(phir[1] - phir[0]);
  double phi_step = phi_span / phi_steps;
  double theta_span = std::abs(thetar[1] - thetar[0]);
  double theta_step = theta_span / theta_steps;
  double access_step = track_length / access_steps;

  return std::visit(
      [&](auto& bField) -> int {
        using field_type =
            typename std::decay_t<decltype(bField)>::element_type;
        if constexpr (!std::is_same_v<field_type, InterpolatedBFieldMap2D> &&
                      !std::is_same_v<field_type, InterpolatedBFieldMap3D>) {
          std::cout << "Bfield map could not be read. Exiting." << std::endl;
          return EXIT_FAILURE;
        } else {
          // Step-wise access pattern
          accessStepWise(*bField, magFieldContext, nEvents, theta_steps,
                         thetar[0], theta_step, phi_steps, phir[0], phi_step,
                         access_steps, access_step);
          // Random access pattern
          accessRandom(*bField, magFieldContext,
                       nEvents * theta_steps * phi_steps * access_steps,
                       track_length);
          return EXIT_SUCCESS;
        }
      },
      bFieldVar);
}

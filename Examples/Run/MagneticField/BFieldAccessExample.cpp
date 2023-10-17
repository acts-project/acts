// This file is part of the Acts project.
//
// Copyright (C) 2017-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/MagneticFieldOptions.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <random>
#include <string>

#include <boost/program_options.hpp>

#if ((BOOST_VERSION / 100) % 1000) <= 71
// Boost <=1.71 and lower do not have progress_display.hpp as a replacement yet
#include <boost/progress.hpp>
using progress_display = boost::progress_display;
#else
// Boost >=1.72 can use this as a replacement
#include <boost/timer/progress_display.hpp>
using progress_display = boost::timer::progress_display;
#endif

/// The main executable
///
/// Creates an InterpolatedBFieldMap from a txt or csv file
/// It then tests random versus stepwise access with the
/// direct getField access and the cell.getField access
/// with cell caching

namespace po = boost::program_options;

using UniformDist = std::uniform_real_distribution<double>;
using RandomEngine = std::mt19937;

void accessStepWise(const Acts::MagneticFieldProvider& bField,
                    const Acts::MagneticFieldContext& bFieldContext,
                    size_t events, size_t theta_steps, double theta_0,
                    double theta_step, size_t phi_steps, double phi_0,
                    double phi_step, size_t access_steps, double access_step) {
  std::cout << "[>>>] Start: step-wise access pattern ... " << std::endl;
  // initialize the field cache
  auto bCache = bField.makeCache(bFieldContext);
  // boost display
  size_t totalSteps = events * theta_steps * phi_steps * access_steps;
  progress_display show_progress(totalSteps);
  // the event loop
  // loop over the events - @todo move to parallel for
  for (size_t ievt = 0; ievt < events; ++ievt) {
    for (size_t itheta = 0; itheta < theta_steps; ++itheta) {
      double theta = theta_0 + itheta * theta_step;
      for (size_t iphi = 0; iphi < phi_steps; ++iphi) {
        double phi = phi_0 + iphi * phi_step;
        // make a direction
        Acts::Vector3 dir(cos(phi) * sin(theta), sin(phi) * sin(theta),
                          cos(theta));
        // check for the current step
        double currentStep = 0.;
        // now step through the magnetic field
        for (size_t istep = 0; istep < access_steps; ++istep) {
          Acts::Vector3 position = currentStep * dir;
          // access the field with the cell
          auto field_from_cache = bField.getField(position, bCache);
          (void)field_from_cache;  // we don't use this explicitly
          // increase the step
          currentStep += access_step;
          // show the progress bar
          ++show_progress;
        }
      }
    }
    std::cout << "[<<<] End result: total steps:" << totalSteps << std::endl;
  }
}

void accessRandom(const Acts::MagneticFieldProvider& bField,
                  const Acts::MagneticFieldContext& bFieldContext,
                  size_t totalSteps, double radius) {
  std::cout << "[>>>] Start: random access pattern ... " << std::endl;
  RandomEngine rng;
  UniformDist xDist(-radius, radius);
  UniformDist yDist(-radius, radius);
  UniformDist zDist(-radius, radius);

  // initialize the field cache
  auto bCache = bField.makeCache(bFieldContext);
  progress_display show_progress(totalSteps);

  // the event loop
  // loop over the events - @todo move to parallel for
  for (size_t istep = 0; istep < totalSteps; ++istep) {
    Acts::Vector3 position(xDist(rng), yDist(rng), zDist(rng));
    // access the field with the cell
    auto field_from_cache = bField.getField(position, bCache);
    (void)field_from_cache;  // we don't use this explicitly
    // show the progress bar
    ++show_progress;
  }
  std::cout << "[<<<] End result: total steps: " << totalSteps << std::endl;
}

/// @brief main executable
///
/// @param argc The argument count
/// @param argv The argument list
int main(int argc, char* argv[]) {
  // Declare the supported program options.
  auto desc = ActsExamples::Options::makeDefaultOptions();
  ActsExamples::Options::addSequencerOptions(desc);
  ActsExamples::Options::addMagneticFieldOptions(desc);
  desc.add_options()(
      "bf-phi-range",
      po::value<ActsExamples::Options::Reals<2>>()->default_value(
          {{-M_PI, M_PI}}),
      "range in which the phi parameter is generated.")(
      "bf-theta-range",
      po::value<ActsExamples::Options::Reals<2>>()->default_value({{0., M_PI}}),
      "range in which the eta parameter is generated.")(
      "bf-phisteps", po::value<size_t>()->default_value(1000),
      "number of steps for the phi parameter.")(
      "bf-thetasteps", po::value<size_t>()->default_value(100),
      "number of steps for the eta parameter.")(
      "bf-accesssteps", po::value<size_t>()->default_value(100),
      "number of steps for magnetic field access.")(
      "bf-tracklength", po::value<double>()->default_value(100.),
      "track length in [mm] magnetic field access.");
  auto vm = ActsExamples::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  // A test magnetic field context
  Acts::MagneticFieldContext magFieldContext = Acts::MagneticFieldContext();

  // TODO
  // Why does this need number-of-events? If it really does emulate
  // per-event access patterns this should be switched to a proper
  // Sequencer-based tool. Otherwise it should be removed.
  auto nEvents = ActsExamples::Options::readSequencerConfig(vm).events.value();
  auto bField = ActsExamples::Options::readMagneticField(vm);

  // Get the phi and eta range
  auto phir = vm["bf-phi-range"].as<ActsExamples::Options::Reals<2>>();
  auto thetar = vm["bf-theta-range"].as<ActsExamples::Options::Reals<2>>();
  // Get the granularity
  size_t phi_steps = vm["bf-phisteps"].as<size_t>();
  size_t theta_steps = vm["bf-thetasteps"].as<size_t>();
  // The defaults
  size_t access_steps = vm["bf-accesssteps"].as<size_t>();
  double track_length =
      vm["bf-tracklength"].as<double>() * Acts::UnitConstants::mm;
  // sort the ranges - and prepare the access grid
  std::sort(phir.begin(), phir.end());
  std::sort(thetar.begin(), thetar.end());
  double phi_span = std::abs(phir[1] - phir[0]);
  double phi_step = phi_span / phi_steps;
  double theta_span = std::abs(thetar[1] - thetar[0]);
  double theta_step = theta_span / theta_steps;
  double access_step = track_length / access_steps;

  accessStepWise(*bField, magFieldContext, nEvents, theta_steps, thetar[0],
                 theta_step, phi_steps, phir[0], phi_step, access_steps,
                 access_step);

  accessRandom(*bField, magFieldContext,
               nEvents * theta_steps * phi_steps * access_steps, track_length);
  return EXIT_SUCCESS;
}

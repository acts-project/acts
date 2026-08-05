// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// Command line options shared by the benchmarks that run on a synthetic event,
/// so that they compare seeders on the same event. Kept out of the library
/// because `boost::program_options` has no business in a shipped Fatras header.

#include "ActsFatras/Synthetic/DetectorLayout.hpp"
#include "ActsFatras/Synthetic/EventGenerator.hpp"

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

namespace ActsTests::SyntheticEventOptions {

/// The options, in the state they have after parsing.
struct Values {
  /// Which shipped detector to generate on
  std::string layout = "itk-pixel";
  /// Number of minimum-bias interactions to overlay
  std::size_t pileup = ActsFatras::Synthetic::GenerationConfig{}.pileup;
  /// Multiplier on both secondary rates the layout is tuned to; a multiplier
  /// and not a value, the two counting in different lengths.
  float secondaryScale = 1.f;
  /// Eta modules per barrel cylinder
  std::uint32_t barrelModules = 1;
  /// Sub-rings each ring the layout declares is split into; one leaves its own
  /// ring structure alone.
  std::uint32_t discRingSplit = 1;
  /// Random seed
  std::uint32_t seed = ActsFatras::Synthetic::EventConfig{}.seed;
};

/// Add the shared options to a description.
/// @param description the description to add them to
/// @param values receives the parsed values; has to outlive the parsing
inline void add(boost::program_options::options_description& description,
                Values& values) {
  namespace po = boost::program_options;
  description.add_options()(
      "layout",
      po::value<std::string>(&values.layout)->default_value(values.layout),
      "detector to generate on: itk-pixel, odd-pixel or generic-pixel")(
      "pileup",
      po::value<std::size_t>(&values.pileup)->default_value(values.pileup),
      "number of overlaid minimum-bias interactions")(
      "secondary-scale",
      po::value<float>(&values.secondaryScale)
          ->default_value(values.secondaryScale),
      "multiplier on the secondary rates the layout is tuned to")(
      "barrel-modules",
      po::value<std::uint32_t>(&values.barrelModules)
          ->default_value(values.barrelModules),
      "eta modules each barrel cylinder is split into")(
      "disc-ring-split",
      po::value<std::uint32_t>(&values.discRingSplit)
          ->default_value(values.discRingSplit),
      "sub-rings each ring the layout declares is split into")(
      "seed",
      po::value<std::uint32_t>(&values.seed)->default_value(values.seed),
      "random seed of the generated event");
}

/// Build the detector the options select.
/// @param values the parsed options
/// @return the layout
inline ActsFatras::Synthetic::DetectorLayout makeLayout(const Values& values) {
  using namespace ActsFatras::Synthetic;

  BarrelEndcapDescription description;
  if (values.layout == "itk-pixel") {
    description = itkPixelDescription();
  } else if (values.layout == "odd-pixel") {
    description = openDataDetectorPixelDescription();
  } else if (values.layout == "generic-pixel") {
    description = genericDetectorPixelDescription();
  } else {
    throw std::invalid_argument("unknown layout '" + values.layout +
                                "', expected itk-pixel, odd-pixel or "
                                "generic-pixel");
  }
  description.barrelModules = values.barrelModules;
  // Subdivide rather than replace: a description carries the real ring
  // structure of its endcap, gaps included, and overriding it with a uniform
  // split would throw that away.
  for (DiscDescription& disc : description.discs) {
    disc.rings = subdivideRings(disc.rings, values.discRingSplit);
  }
  return ActsFatras::Synthetic::makeLayout(description);
}

/// Build the event configuration the options select. It follows the layout:
/// running the ODD layout with the ITk numbers would generate a quarter too
/// many space points.
///
/// @param values the parsed options
/// @return the configuration
inline ActsFatras::Synthetic::EventConfig makeConfig(const Values& values) {
  using namespace ActsFatras::Synthetic;

  EventConfig config;
  if (values.layout == "itk-pixel") {
    config = EventConfig::itkPixelTtbarPu200();
  } else {
    // The Generic detector has no full simulation of its own to be fitted
    // against, and the ODD is its descendant and close enough in size.
    config = EventConfig::openDataDetectorTtbarPu200();
  }

  config.generation.pileup = values.pileup;
  config.seed = values.seed;
  SecondaryConfig& secondaries = config.simulation.secondaries;
  secondaries.electronRate *= values.secondaryScale;
  secondaries.nuclearRate *= values.secondaryScale;
  return config;
}

}  // namespace ActsTests::SyntheticEventOptions

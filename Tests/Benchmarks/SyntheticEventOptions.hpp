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

#include "ActsFatras/Json/DataDirectory.hpp"
#include "ActsFatras/Json/DetectorDescriptionJsonConverter.hpp"
#include "ActsFatras/Json/EventConfigJsonConverter.hpp"
#include "ActsFatras/Synthetic/DetectorLayout.hpp"
#include "ActsFatras/Synthetic/EventConfig.hpp"
#include "ActsFatras/Synthetic/EventGenerator.hpp"

#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

namespace ActsTests::SyntheticEventOptions {

/// The options, in the state they have after parsing.
struct Values {
  /// Which detector to generate on: one that ships in `Fatras/data`, or the
  /// path that one's files share; see `layoutPath`
  std::string layout = "itk";
  /// Which of its subsystems to build, comma separated; empty builds all of
  /// them
  std::string subsystems;
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
      "detector to generate on: itk or odd, or the path its files share, e.g. "
      "/tmp/mine for /tmp/mine-description.json")(
      "subsystems", po::value<std::string>(&values.subsystems),
      "subsystems of it to build, comma separated; all of them by default")(
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

/// One of a detector's files. A detector is three of them that share a prefix,
/// so `--layout` names that prefix: `itk` for the shipped ones, or a path like
/// `/tmp/mine` for `/tmp/mine-description.json` beside its material and its
/// configuration.
///
/// @param values the parsed options
/// @param suffix which of its files, e.g. `-description.json`
/// @return the path to read
inline std::filesystem::path layoutPath(const Values& values,
                                        const std::string& suffix) {
  const std::filesystem::path given{values.layout};
  if (given.has_parent_path()) {
    return given.parent_path() / (given.filename().string() + suffix);
  }
  return ActsFatras::dataPath(values.layout + suffix);
}

/// Build the detector the options select, out of the description and material
/// it ships as.
///
/// @param values the parsed options
/// @return the layout
/// @throws std::runtime_error if there is no such detector
inline ActsFatras::Synthetic::DetectorLayout makeLayout(const Values& values) {
  using namespace ActsFatras::Synthetic;

  DetectorDescription description =
      readDetectorDescription(layoutPath(values, "-description.json"));
  decorate(description,
           readMaterialDecoration(layoutPath(values, "-material.json")));

  if (!values.subsystems.empty()) {
    std::vector<std::string> names;
    boost::split(names, values.subsystems, boost::is_any_of(","));
    description = selectSubsystems(description, names);
  }

  for (SubsystemDescription& subsystem : description.subsystems) {
    for (BarrelDescription& barrel : subsystem.barrels) {
      for (CylinderDescription& cylinder : barrel.cylinders) {
        cylinder.modules = values.barrelModules;
      }
    }
    for (EndcapDescription& endcap : subsystem.endcaps) {
      for (DiscDescription& disc : endcap.discs) {
        // Subdivide rather than replace: a description carries the real ring
        // structure of its endcap, gaps included, and overriding it with a
        // uniform split would throw that away.
        disc.rings = subdivideRings(disc.rings, values.discRingSplit);
      }
    }
  }

  return ActsFatras::Synthetic::makeLayout(description);
}

/// Build the event configuration fitted to that detector. It follows the
/// layout: running the ODD with the ITk numbers would generate a quarter too
/// many space points.
///
/// @param values the parsed options
/// @return the configuration
inline ActsFatras::Synthetic::EventConfig makeConfig(const Values& values) {
  using namespace ActsFatras::Synthetic;

  EventConfig config = readEventConfig(layoutPath(values, "-ttbar-pu200.json"));
  config.generation.pileup = values.pileup;
  config.seed = values.seed;
  SecondaryConfig& secondaries = config.simulation.secondaries;
  secondaries.electronRate *= values.secondaryScale;
  secondaries.nuclearRate *= values.secondaryScale;
  return config;
}

}  // namespace ActsTests::SyntheticEventOptions

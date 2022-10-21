// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Options/TruthSeedSelectorOptions.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <boost/program_options.hpp>

namespace ActsExamples {
namespace Options {

void addTruthSeedSelectorOptions(Options::Description& desc) {
  using boost::program_options::value;
  using Options::Interval;

  auto opt = desc.add_options();
  opt("select-rho-mm", value<Interval>()->value_name("MIN:MAX"),
      "Select particle transverse distance to the origin in mm");
  opt("select-z-mm", value<Interval>()->value_name("MIN:MAX"),
      "Select particle longitudinal distance range to the origin in mm");
  opt("select-phi-degree", value<Interval>()->value_name("MIN:MAX"),
      "Select particle direction angle in the transverse plane in degree");
  opt("select-eta", value<Interval>()->value_name("MIN:MAX"),
      "Select particle pseudo-rapidity");
  opt("select-abseta", value<Interval>()->value_name("MIN:MAX"),
      "Select particle absolute pseudo-rapidity");
  opt("select-pt-gev", value<Interval>()->value_name("MIN:MAX"),
      "Select particle transverse momentum in GeV");
  opt("select-min-hits", value<size_t>()->default_value(3),
      "Select particle minimum hits");
}

ActsExamples::TruthSeedSelector::Config readTruthSeedSelectorConfig(
    const Options::Variables& vars) {
  using namespace Acts::UnitLiterals;

  // Set boundary values if the given config exists
  auto extractInterval = [&](const char* name, auto unit, auto& lower,
                             auto& upper) {
    if (vars[name].empty()) {
      return;
    }
    auto interval = vars[name].as<Options::Interval>();
    lower = interval.lower.value_or(lower) * unit;
    upper = interval.upper.value_or(upper) * unit;
  };

  ActsExamples::TruthSeedSelector::Config cfg;
  extractInterval("select-rho-mm", 1_mm, cfg.rhoMin, cfg.rhoMax);
  extractInterval("select-z-mm", 1_mm, cfg.zMin, cfg.zMax);
  extractInterval("select-phi-degree", 1_degree, cfg.phiMin, cfg.phiMax);
  extractInterval("select-eta", 1.0, cfg.etaMin, cfg.etaMax);
  extractInterval("select-abseta", 1.0, cfg.absEtaMin, cfg.absEtaMax);
  extractInterval("select-pt-gev", 1_GeV, cfg.ptMin, cfg.ptMax);
  cfg.nHitsMin = vars["select-min-hits"].as<size_t>();
  return cfg;
}

}  // namespace Options
}  // namespace ActsExamples

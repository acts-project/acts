// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Options/ParticleGunOptions.hpp"

#include "Acts/Utilities/Units.hpp"
#include "ActsExamples/Generators/MultiplicityGenerators.hpp"
#include "ActsExamples/Generators/ParametricParticleGenerator.hpp"
#include "ActsExamples/Generators/VertexGenerators.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <boost/program_options.hpp>

void ActsExamples::Options::addParticleGunOptions(Description& desc) {
  using boost::program_options::bool_switch;
  using boost::program_options::value;

  auto opt = desc.add_options();
  opt("pg-vertex-xy-std-mm", value<double>()->default_value(0.0),
      "Transverse vertex standard deviation in mm");
  opt("pg-vertex-z-std-mm", value<double>()->default_value(0.0),
      "Longitudinal vertex standard deviation in mm");
  opt("pg-vertex-t-std-ns", value<double>()->default_value(0.0),
      "Temporal vertex standard deviation in ns");
  opt("pg-phi-degree",
      value<Interval>()->value_name("MIN:MAX")->default_value({0.0, 360.0}),
      "Transverse direction angle generation range in degree");
  opt("pg-eta",
      value<Interval>()->value_name("MIN:MAX")->default_value({-4.0, 4.0}),
      "Pseudo-rapidity generation range");
  opt("pg-p-gev",
      value<Interval>()->value_name("MIN:MAX")->default_value({1.0, 10.0}),
      "Absolute momentum generation range in GeV");
  opt("pg-pdg", value<int32_t>()->default_value(Acts::PdgParticle::eMuon),
      "PDG number of the particle, will be adjusted for charge flip.");
  opt("pg-randomize-charge", bool_switch(),
      "Flip the charge and change the PDG number accordingly.");
  opt("pg-nparticles", value<size_t>()->default_value(1u),
      "Number of generated particles");
}

ActsExamples::EventGenerator::Config
ActsExamples::Options::readParticleGunOptions(const Variables& vars) {
  using namespace Acts::UnitLiterals;

  // access user config w/ unit conversion
  auto getValue = [&](const char* name, auto unit) {
    return vars[name].as<double>() * unit;
  };
  auto getRange = [&](const char* name, auto unit, auto& lower, auto& upper) {
    auto interval = vars[name].as<Options::Interval>();
    lower = interval.lower.value() * unit;
    upper = interval.upper.value() * unit;
  };

  GaussianVertexGenerator vertexGen;
  vertexGen.stddev[Acts::ePos0] = getValue("pg-vertex-xy-std-mm", 1_mm);
  vertexGen.stddev[Acts::ePos1] = getValue("pg-vertex-xy-std-mm", 1_mm);
  vertexGen.stddev[Acts::ePos2] = getValue("pg-vertex-z-std-mm", 1_mm);
  vertexGen.stddev[Acts::eTime] = getValue("pg-vertex-t-std-ns", 1_ns);

  ParametricParticleGenerator::Config pgCfg;
  getRange("pg-phi-degree", 1_degree, pgCfg.phiMin, pgCfg.phiMax);
  // user config sets eta but the generator takes theta
  double etaMin, etaMax;
  getRange("pg-eta", 1.0, etaMin, etaMax);
  pgCfg.thetaMin = 2 * std::atan(std::exp(-etaMin));
  pgCfg.thetaMax = 2 * std::atan(std::exp(-etaMax));
  getRange("pg-p-gev", 1_GeV, pgCfg.pMin, pgCfg.pMax);
  pgCfg.pdg =
      static_cast<Acts::PdgParticle>(vars["pg-pdg"].template as<int32_t>());
  pgCfg.randomizeCharge = vars["pg-randomize-charge"].template as<bool>();
  pgCfg.numParticles = vars["pg-nparticles"].as<size_t>();

  EventGenerator::Config cfg;
  cfg.generators = {
      {FixedMultiplicityGenerator{1}, std::move(vertexGen),
       ParametricParticleGenerator(pgCfg)},
  };

  return cfg;
}

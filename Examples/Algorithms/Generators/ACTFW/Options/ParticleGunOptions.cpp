// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Options/ParticleGunOptions.hpp"

#include "ACTFW/Generators/MultiplicityGenerators.hpp"
#include "ACTFW/Generators/ParametricProcessGenerator.hpp"
#include "ACTFW/Generators/VertexGenerators.hpp"
#include "ACTFW/Utilities/Options.hpp"
#include "Acts/Utilities/Units.hpp"

void FW::Options::addParticleGunOptions(
    boost::program_options::options_description& opt) {
  using namespace boost::program_options;

  opt.add_options()("pg-nparticles", value<size_t>()->default_value(1.),
                    "number of particles.")(
      "pg-d0-range", value<read_range>()->multitoken()->default_value({0., 0.}),
      "range in which the d0 parameter is simulated in [mm]. Please hand"
      "over by simply seperating the values by space")(
      "pg-z0-range", value<read_range>()->multitoken()->default_value({0., 0.}),
      "range in which the z0 parameter is simulated in [mm]. Please hand"
      "over by simply seperating the values by space")(
      "pg-t0-range", value<read_range>()->multitoken()->default_value({0., 0.}),
      "range in which the t0 parameter is simulated in [ns]. Please hand"
      "over by simply seperating the values by space")(
      "pg-phi-range",
      value<read_range>()->multitoken()->default_value({-M_PI, M_PI}),
      "range in which the phi0 parameter is simulated. Please hand over by "
      "simply seperating the values by space")(
      "pg-eta-range",
      value<read_range>()->multitoken()->default_value({-4., 4.}),
      "range in which the eta parameter is simulated. Please hand over by "
      "simply seperating the values by space")(
      "pg-pt-range", value<read_range>()->multitoken()->default_value({1, 100}),
      "range in which the pt in [GeV] parameter is simulated. Please hand "
      "over by simply seperating the values by space")(
      "pg-pdg", value<int32_t>()->default_value(Acts::PdgParticle::eMuon),
      "PDG number of the particle, will be adjusted for charge flip.")(
      "pg-randomize-charge", bool_switch(),
      "flip the charge (and change PDG accordingly).");
}

FW::EventGenerator::Config FW::Options::readParticleGunOptions(
    const boost::program_options::variables_map& vm) {
  using namespace Acts::UnitLiterals;

  // read the range as vector (missing istream for std::array)
  auto d0 = vm["pg-d0-range"].template as<read_range>();
  auto z0 = vm["pg-z0-range"].template as<read_range>();
  auto t0 = vm["pg-t0-range"].template as<read_range>();
  auto phi = vm["pg-phi-range"].template as<read_range>();
  auto eta = vm["pg-eta-range"].template as<read_range>();
  auto pt = vm["pg-pt-range"].template as<read_range>();

  ParametricProcessGenerator::Config pgCfg;
  pgCfg.numParticles = vm["pg-nparticles"].template as<size_t>();
  pgCfg.d0Range = {{d0[0] * 1_mm, d0[1] * 1_mm}};
  pgCfg.z0Range = {{z0[0] * 1_mm, z0[1] * 1_mm}};
  pgCfg.t0Range = {{t0[0] * 1_ns, t0[1] * 1_ns}};
  pgCfg.phiRange = {{phi[0], phi[1]}};
  pgCfg.etaRange = {{eta[0], eta[1]}};
  pgCfg.ptRange = {{pt[0] * 1_GeV, pt[1] * 1_GeV}};
  pgCfg.pdg =
      static_cast<Acts::PdgParticle>(vm["pg-pdg"].template as<int32_t>());
  pgCfg.randomizeCharge = vm["pg-randomize-charge"].template as<bool>();

  EventGenerator::Config cfg;
  cfg.generators = {
      {FixedMultiplicityGenerator{1},
       FixedVertexGenerator{{0.0, 0.0, 0.0, 0.0}},
       ParametricProcessGenerator{pgCfg}},
  };

  return cfg;
}

// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Options/Pythia8Options.hpp"

#include <boost/program_options.hpp>

#include "ACTFW/Generators/MultiplicityGenerators.hpp"
#include "ACTFW/Generators/Pythia8ProcessGenerator.hpp"
#include "ACTFW/Generators/VertexGenerators.hpp"
#include "Acts/Utilities/PdgParticle.hpp"
#include "Acts/Utilities/Units.hpp"

void FW::Options::addPythia8Options(
    boost::program_options::options_description& opt) {
  using namespace boost::program_options;

  opt.add_options()("evg-cms-energy", value<double>()->default_value(14000.),
                    "Center-of-mass energy collision in GeV")(
      "evg-beam0", value<int32_t>()->default_value(Acts::PdgParticle::eProton),
      "PDG number of the first beam particle")(
      "evg-beam1", value<int32_t>()->default_value(Acts::PdgParticle::eProton),
      "PDG number of the second beam particle")(
      "evg-hard-process",
      value<std::string>()->default_value("HardQCD:all = on"),
      "Pythia8 process string for the hard scatter")(
      "evg-pileup-process",
      value<std::string>()->default_value("SoftQCD:all = on"),
      "Pythi8 process string for the pile-up")(
      "evg-pileup", value<size_t>()->default_value(200),
      "Number of instantaneous pile-up events")(
      "evg-vertex-xy-std", value<double>()->default_value(0.015),
      "Transverse vertex standard deviation in mm")(
      "evg-vertex-z-std", value<double>()->default_value(55.5),
      "Longitudinal vertex standard deviation in mm")(
      "evg-vertex-t-std", value<double>()->default_value(0.08),
      "Temporal vertex standard deviation in ns")(
      "evg-shuffle", bool_switch(), "Randomnly shuffle the vertex order.");
}

FW::EventGenerator::Config FW::Options::readPythia8Options(
    const boost::program_options::variables_map& vm, Acts::Logging::Level lvl) {
  Pythia8Generator::Config hard;
  hard.pdgBeam0 = static_cast<Acts::PdgParticle>(vm["evg-beam0"].as<int32_t>());
  hard.pdgBeam1 = static_cast<Acts::PdgParticle>(vm["evg-beam1"].as<int32_t>());
  hard.cmsEnergy = vm["evg-cms-energy"].as<double>();
  hard.settings = {vm["evg-hard-process"].as<std::string>()};

  Pythia8Generator::Config pileup;
  pileup.pdgBeam0 =
      static_cast<Acts::PdgParticle>(vm["evg-beam0"].as<int32_t>());
  pileup.pdgBeam1 =
      static_cast<Acts::PdgParticle>(vm["evg-beam1"].as<int32_t>());
  pileup.cmsEnergy = vm["evg-cms-energy"].as<double>();
  pileup.settings = {vm["evg-pileup-process"].as<std::string>()};

  auto mu = vm["evg-pileup"].as<size_t>();
  auto vtxStdXY =
      vm["evg-vertex-xy-std"].as<double>() * Acts::UnitConstants::mm;
  auto vtxStdZ = vm["evg-vertex-z-std"].as<double>() * Acts::UnitConstants::mm;
  auto vtxStdT = vm["evg-vertex-t-std"].as<double>() * Acts::UnitConstants::ns;

  EventGenerator::Config cfg;
  cfg.generators = {
      {FixedMultiplicityGenerator{1},
       GaussianVertexGenerator{{vtxStdXY, vtxStdXY, vtxStdZ, vtxStdT}},
       Pythia8Generator::makeFunction(hard, lvl)},
      {PoissonMultiplicityGenerator{mu},
       GaussianVertexGenerator{{vtxStdXY, vtxStdXY, vtxStdZ, vtxStdT}},
       Pythia8Generator::makeFunction(hard, lvl)},
  };
  cfg.shuffle = vm["evg-shuffle"].as<bool>();

  return cfg;
}

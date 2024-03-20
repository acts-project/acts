// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Options/ParticleSmearingOptions.hpp"

void ActsExamples::Options::addParticleSmearingOptions(Description& desc) {
  using boost::program_options::value;
  using Options::Reals;

  auto opt = desc.add_options();
  opt("smear-sigma-D0", value<Reals<3>>()->default_value({{20, 30, 0.3}}),
      "Smear the initial Pt-dependent d0 in perigee frame with a_0[um] + "
      "a_1[um]*exp(-1.*abs(a_2[1/GeV])*pt)");
  opt("smear-sigma-Z0", value<Reals<3>>()->default_value({{20, 30, 0.3}}),
      "Smear the initial Pt-dependent z0 in perigee frame with a_0[um] + "
      "a_1[um]*exp(-1.*abs(a_2[1/GeV])*pt)");
  opt("smear-sigma-T0", value<double>()->default_value(1),
      "Smear the initial time in ns");
  opt("smear-sigma-momentum", value<Reals<3>>()->default_value({{1, 1, 0.05}}),
      "Smear the initial phi (degree), theta (degree) and momentum (relative)");
  opt("smear-initial-variance-inflation",
      value<Reals<6>>()->default_value({{1., 1., 1., 1., 1., 1.}}),
      "Inflate the initial covariance matrix");
}

ActsExamples::ParticleSmearing::Config
ActsExamples::Options::readParticleSmearingOptions(
    const ActsExamples::Options::Variables& vars) {
  using namespace ActsExamples;
  using namespace Acts::UnitConstants;
  using Options::Reals;

  ParticleSmearing::Config cfg;

  auto sigmaD0Opts = vars["smear-sigma-D0"].template as<Reals<3>>();
  auto sigmaZ0Opts = vars["smear-sigma-Z0"].template as<Reals<3>>();
  auto sigmaMomOpts = vars["smear-sigma-momentum"].template as<Reals<3>>();
  cfg.sigmaD0 = sigmaD0Opts[0] * Acts::UnitConstants::um;
  cfg.sigmaD0PtA = sigmaD0Opts[1] * Acts::UnitConstants::um;
  cfg.sigmaD0PtB = sigmaD0Opts[2] / Acts::UnitConstants::GeV;
  cfg.sigmaZ0 = sigmaZ0Opts[0] * Acts::UnitConstants::um;
  cfg.sigmaZ0PtA = sigmaZ0Opts[1] * Acts::UnitConstants::um;
  cfg.sigmaZ0PtB = sigmaZ0Opts[2] / Acts::UnitConstants::GeV;
  cfg.sigmaT0 = vars["smear-sigma-T0"].as<double>() * Acts::UnitConstants::ns;
  cfg.sigmaPhi = sigmaMomOpts[0] * Acts::UnitConstants::degree;
  cfg.sigmaTheta = sigmaMomOpts[1] * Acts::UnitConstants::degree;
  cfg.sigmaPRel = sigmaMomOpts[2];
  auto varInflation =
      vars["smear-initial-variance-inflation"].template as<Reals<6>>();
  cfg.initialVarInflation = {varInflation[0], varInflation[1], varInflation[2],
                             varInflation[3], varInflation[4], varInflation[5]};
  return cfg;
}

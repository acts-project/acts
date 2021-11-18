// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/ParticleSmearingOptions.hpp"

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
  opt("smear-sigma-momentum", value<Reals<3>>()->default_value({{1, 1, 0.1}}),
      "Smear the initial phi (degree), theta (degree) and momentum (relative)");
}
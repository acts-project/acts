// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Fatras/FatrasOptions.hpp"

#include <string>

void FW::Options::addFatrasOptions(FW::Options::Description& desc) {
  using boost::program_options::bool_switch;
  using boost::program_options::value;

  auto opt = desc.add_options();
  opt("fatras-pmin-gev", value<double>()->default_value(0.5),
      "Minimum momentum for simulated particles in GeV");
  opt("fatras-em-scattering", value<bool>()->default_value(true),
      "Simulate multiple scattering of charged particles");
  opt("fatras-em-ionisation", value<bool>()->default_value(true),
      "Simulate ionisiation/excitation energy loss of charged particles");
  opt("fatras-em-radiation", value<bool>()->default_value(true),
      "Simulate radiative energy loss of charged particles");
  opt("fatras-hits",
      value<std::string>()
          ->value_name("none|sensitive|material|all")
          ->default_value("sensitive"),
      "Which surfaces should record charged particle hits");
}

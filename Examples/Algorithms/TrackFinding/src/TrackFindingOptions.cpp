// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/TrackFinding/TrackFindingOptions.hpp"

#include <boost/program_options.hpp>
#include <string>

void FW::Options::addTrackFindingOptions(FW::Options::Description& desc) {
  using boost::program_options::value;

  auto opt = desc.add_options();
  opt("ckf-slselection-chi2max", value<double>()->default_value(15),
      "Global criteria of maximum chi2 for CKF source link selection");
  opt("ckf-slselection-nmax", value<int>()->default_value(10),
      "Global criteria of maximum number of source link candidates on a "
      "surface for CKF source link selection");
}

FW::TrackFindingAlgorithm::Config FW::Options::readTrackFindingConfig(
    const FW::Options::Variables& variables) {
  using Config = typename FW::TrackFindingAlgorithm::Config;
  Config tfAlgCfg;

  tfAlgCfg.sourcelinkSelectorCfg.globalChi2CutOff =
      variables["ckf-slselection-chi2max"].template as<double>();
  tfAlgCfg.sourcelinkSelectorCfg.globalNumSourcelinksCutOff =
      variables["ckf-slselection-nmax"].template as<int>();

  return tfAlgCfg;
}

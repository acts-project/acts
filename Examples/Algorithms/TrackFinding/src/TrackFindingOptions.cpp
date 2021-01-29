// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/TrackFindingOptions.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <string>

#include <boost/program_options.hpp>

void ActsExamples::Options::addTrackFindingOptions(
    ActsExamples::Options::Description& desc) {
  using boost::program_options::value;

  auto opt = desc.add_options();
  opt("ckf-selection-chi2max", value<double>()->default_value(15),
      "Global criteria of maximum chi2 for CKF measurement selection");
  opt("ckf-selection-nmax", value<size_t>()->default_value(10),
      "Global criteria of maximum number of measurement candidates on a "
      "surface for CKF measurement selection");
}

ActsExamples::TrackFindingAlgorithm::Config
ActsExamples::Options::readTrackFindingConfig(
    const ActsExamples::Options::Variables& variables) {
  auto chi2Max = variables["ckf-selection-chi2max"].template as<double>();
  auto nMax = variables["ckf-selection-nmax"].template as<size_t>();

  // config is a GeometryHierarchyMap with just the global default
  TrackFindingAlgorithm::Config cfg;
  cfg.measurementSelectorCfg = {
      {Acts::GeometryIdentifier(), {chi2Max, nMax}},
  };
  return cfg;
}

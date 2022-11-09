// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Options/TrackFindingOptions.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <string>

#include <boost/program_options.hpp>

void ActsExamples::Options::addTrackFindingOptions(
    ActsExamples::Options::Description& desc) {
  using boost::program_options::value;

  auto opt = desc.add_options();
  opt("ckf-selection-abseta-bins", value<VariableReals>()->default_value({{}}),
      "bins in |eta| to specify variable selections");
  opt("ckf-selection-chi2max", value<VariableReals>()->default_value({{15}}),
      "Maximum chi2 for CKF measurement selection "
      "(specify multiple values with --ckf-selection-abseta-bins)");
  opt("ckf-selection-nmax", value<VariableIntegers>()->default_value({{10}}),
      "Maximum number of measurement candidates on a "
      "surface for CKF measurement selection "
      "(specify multiple values with --ckf-selection-abseta-bins)");
  opt("ckf-initial-variance-inflation",
      value<Reals<6>>()->default_value({{1., 1., 1., 1., 1., 1.}}),
      "Inflation factor for the initial variances in the CKF search, must be "
      "of form i:j:k:l:m:n.");
}

ActsExamples::TrackFindingAlgorithm::Config
ActsExamples::Options::readTrackFindingConfig(
    const ActsExamples::Options::Variables& variables) {
  auto etaBins = variables["ckf-selection-abseta-bins"]
                     .template as<VariableReals>()
                     .values;
  auto chi2Max =
      variables["ckf-selection-chi2max"].template as<VariableReals>().values;
  auto nMax =
      variables["ckf-selection-nmax"].template as<VariableIntegers>().values;

  // config is a GeometryHierarchyMap with just the global default
  TrackFindingAlgorithm::Config cfg;
  cfg.measurementSelectorCfg = {
      {Acts::GeometryIdentifier(),
       {etaBins, chi2Max, {nMax.begin(), nMax.end()}}},
  };
  return cfg;
}

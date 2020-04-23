// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/program_options.hpp>
#include <cstdlib>
#include <iostream>
#include <utility>

#include "ACTFW/Utilities/Options.hpp"
#include "Acts/Plugins/TGeo/TGeoLayerBuilder.hpp"
#include "Acts/Utilities/Units.hpp"

namespace po = boost::program_options;

namespace au = Acts::units;

namespace FW {

namespace Options {

/// The root detecotr options, the are prefixes with geo-tgeo
///
/// @tparam options_t Type of the options object (bound to boost API)
///
/// @param opt The provided object, where root specific options are attached
template <typename options_t>
void addTGeoGeometryOptions(options_t& opt) {
  opt.add_options()("geo-tgeo-filename",
                    po::value<std::string>()->default_value(""),
                    "Root file name.")(
      "geo-tgeo-worldvolume", po::value<std::string>()->default_value(""),
      "Root world volume to start search from.")(
      "geo-tgeo-unitScalor", po::value<double>()->default_value(10.),
      "Unit scalor from ROOT to Acts.")(
      "geo-tgeo-nlayers",
      po::value<read_series>()->multitoken()->default_value({}),
      "Number of layers on the negative side.")(
      "geo-tgeo-clayers",
      po::value<read_series>()->multitoken()->default_value({}),
      "Number of layers in the barrel.")(
      "geo-tgeo-players",
      po::value<read_series>()->multitoken()->default_value({}),
      "Number of layers on the positive side.")(
      "geo-tgeo-nlayernames",
      po::value<read_strings>()->multitoken()->default_value({}),
      "Name identifier for negative layer objects, odered along the series.")(
      "geo-tgeo-clayernames",
      po::value<read_strings>()->multitoken()->default_value({}),
      "Name identifier for central layer objects, odered along the series.")(
      "geo-tgeo-playernames",
      po::value<read_strings>()->multitoken()->default_value({}),
      "Name identifier for positive layer objects, odered along the series.")(
      "geo-tgeo-nmodulenames",
      po::value<read_strings>()->multitoken()->default_value({}),
      "Name identifier for negative sensitive objects, odered along the "
      "series.")("geo-tgeo-cmodulenames",
                 po::value<read_strings>()->multitoken()->default_value({}),
                 "Name identifier for central sensitive objects, odered "
                 "along the series.")(
      "geo-tgeo-pmodulenames",
      po::value<read_strings>()->multitoken()->default_value({}),
      "Name identifier for positive sensitive objects, odered along the "
      "series.")("geo-tgeo-nmoduleaxes",
                 po::value<read_strings>()->multitoken()->default_value({}),
                 "Axes definition for negative sensitive objects, odered "
                 "along the series.")(
      "geo-tgeo-cmoduleaxes",
      po::value<read_strings>()->multitoken()->default_value({}),
      "Axes definition for central sensitive objects, odered along the "
      "series.")("geo-tgeo-pmoduleaxes",
                 po::value<read_strings>()->multitoken()->default_value({}),
                 "Axes definition for positive sensitive objects, odered "
                 "along the series.");
}

/// @brief Read the specific options for the ROOT detector
///
/// @tparam variable_map_t Type of the variable matp to read out
///
/// @param vm The variable map containing the options
/// @param layerCreator The provide layer creator
///
/// @return a configuration object for a TGeoLayerBuilder
template <typename variable_map_t>
std::vector<Acts::TGeoLayerBuilder::Config> readTGeoLayerBuilderConfigs(
    const variable_map_t& vm,
    std::shared_ptr<const Acts::LayerCreator> layerCreator) {
  std::vector<Acts::TGeoLayerBuilder::Config> detLayerConfigs;

  // general stuff
  read_strings subdetectors =
      vm["geo-subdetectors"].template as<read_strings>();
  double unitScalor = vm["geo-tgeo-unitScalor"].template as<double>();
  // these define a series, such as 1, 3, 4
  read_series nlayers = vm["geo-tgeo-nlayers"].template as<read_series>();
  read_series clayers = vm["geo-tgeo-clayers"].template as<read_series>();
  read_series players = vm["geo-tgeo-players"].template as<read_series>();
  // these are going continously through the series, such as i in [ 1 + 3 + 4
  // ]
  read_strings nlayernames =
      vm["geo-tgeo-nlayernames"].template as<read_strings>();
  read_strings clayernames =
      vm["geo-tgeo-clayernames"].template as<read_strings>();
  read_strings playernames =
      vm["geo-tgeo-playernames"].template as<read_strings>();
  read_strings nsensitivenames =
      vm["geo-tgeo-nmodulenames"].template as<read_strings>();
  read_strings csensitivenames =
      vm["geo-tgeo-cmodulenames"].template as<read_strings>();
  read_strings psensitivenames =
      vm["geo-tgeo-pmodulenames"].template as<read_strings>();
  read_strings nsensitiveaxes =
      vm["geo-tgeo-nmoduleaxes"].template as<read_strings>();
  read_strings csensitiveaxes =
      vm["geo-tgeo-cmoduleaxes"].template as<read_strings>();
  read_strings psensitiveaxes =
      vm["geo-tgeo-pmoduleaxes"].template as<read_strings>();
  // todo consistency checks

  // total coutners
  int tin = 0;
  int tic = 0;
  int tip = 0;
  // prepare the TGeoLayerBuilder::Configs
  for (size_t idet = 0; idet < clayers.size(); ++idet) {
    Acts::TGeoLayerBuilder::Config layerBuilderConfig;
    // get the current numbers
    int cn = nlayers[idet];
    int cc = clayers[idet];
    int cp = players[idet];
    // fill the configs - for negative
    for (int in = 0; in < cn; ++in, ++tin) {
      Acts::TGeoLayerBuilder::LayerConfig lConfig;
      lConfig.layerName = nlayernames[tin];
      lConfig.sensorName = nsensitivenames[tin];
      lConfig.localAxes = nsensitiveaxes[tin];
      layerBuilderConfig.layerConfigurations[0].push_back(lConfig);
    }
    // fill the configs - for barrel
    for (int ic = 0; ic < cc; ++ic, ++tic) {
      Acts::TGeoLayerBuilder::LayerConfig lConfig;
      lConfig.layerName = clayernames[tic];
      lConfig.sensorName = csensitivenames[tic];
      lConfig.localAxes = csensitiveaxes[tic];
      layerBuilderConfig.layerConfigurations[1].push_back(lConfig);
    }
    // fill the configs - for positive
    for (int ip = 0; ip < cp; ++ip, ++tip) {
      Acts::TGeoLayerBuilder::LayerConfig lConfig;
      lConfig.layerName = playernames[tip];
      lConfig.sensorName = psensitivenames[tip];
      lConfig.localAxes = psensitiveaxes[tip];
      layerBuilderConfig.layerConfigurations[2].push_back(lConfig);
    }
    // set the scale and the layer creator
    layerBuilderConfig.configurationName = subdetectors[idet];
    layerBuilderConfig.unit = unitScalor;
    layerBuilderConfig.layerCreator = layerCreator;
    // now add it to the configs
    detLayerConfigs.push_back(layerBuilderConfig);
  }
  return detLayerConfigs;
}
}  // namespace Options
}  // namespace FW

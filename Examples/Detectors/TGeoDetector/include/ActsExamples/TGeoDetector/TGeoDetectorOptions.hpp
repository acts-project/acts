// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/TGeo/TGeoLayerBuilder.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <utility>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

namespace ActsExamples {

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
      "geo-tgeo-unit-scalor", po::value<double>()->default_value(10.),
      "Unit scalor from ROOT to Acts.")(
      "geo-tgeo-bp-parameters",
      po::value<read_range>()->multitoken()->default_value({}),
      "Potential beam pipe parameters {r, z, t} in [mm].")(
      "geo-tgeo-nlayers", po::value<read_series>()->default_value({}),
      "Number of layers on the negative side.")(
      "geo-tgeo-clayers", po::value<read_series>()->default_value({}),
      "Number of layers in the barrel.")(
      "geo-tgeo-players", po::value<read_series>()->default_value({}),
      "Number of layers on the positive side.")(
      "geo-tgeo-ringlayout", po::value<read_series>()->default_value({}),
      "Indicator if ring layout is present.")(
      "geo-tgeo-ringtolerance", po::value<read_range>()->default_value({}),
      "Tolerance for ring layout detection in [mm].")(
      "geo-tgeo-nvolume-name", po::value<read_strings>()->default_value({}),
      "Name identifier of the volume for searching negative layers.")(
      "geo-tgeo-cvolume-name", po::value<read_strings>()->default_value({}),
      "Name identifier of the volume for searching central layers.")(
      "geo-tgeo-pvolume-name", po::value<read_strings>()->default_value({}),
      "Name identifier of the volume for searching positive layers.")(
      "geo-tgeo-nmodule-name", po::value<read_strings>()->default_value({}),
      "Name identifier for negative sensitive objects, odered along the "
      "series.")("geo-tgeo-cmodule-name",
                 po::value<read_strings>()->default_value({}),
                 "Name identifier for central sensitive objects, odered "
                 "along the series.")(
      "geo-tgeo-pmodule-name", po::value<read_strings>()->default_value({}),
      "Name identifier for positive sensitive objects, odered along the "
      "series.")("geo-tgeo-nlayer-r-range",
                 po::value<std::vector<Interval>>()->default_value({}),
                 "Radial range(s) for negative layers "
                 "to restrict the module parsing (optional).")(
      "geo-tgeo-clayer-r-range",
      po::value<std::vector<Interval>>()->default_value({}),
      "Radial range(s) for central layers "
      "to restrict the module parsing (optional).")(
      "geo-tgeo-player-r-range",
      po::value<std::vector<Interval>>()->default_value({}),
      "Radial range(s) for positive layers "
      "to restrict the module parsing (optional).")(
      "geo-tgeo-nlayer-z-range",
      po::value<std::vector<Interval>>()->default_value({}),
      "Longitudinal range(s) for negative layers "
      "to restrict the module parsing (optional).")(
      "geo-tgeo-clayer-z-range",
      po::value<std::vector<Interval>>()->default_value({}),
      "Longitudinal range(s) for central layers "
      "to restrict the module parsing (optional).")(
      "geo-tgeo-player-z-range",
      po::value<std::vector<Interval>>()->default_value({}),
      "Longitudinal range(s) for positive layers "
      "to restrict the module parsing (optional).")(
      "geo-tgeo-nlayer-r-split", po::value<read_range>()->default_value({}),
      "R-tolerances (if > 0.) that triggers splitting "
      " of collected surfaces into different negative layers.")(
      "geo-tgeo-nlayer-z-split", po::value<read_range>()->default_value({}),
      "Z-tolerances (if > 0.) that triggers splitting "
      " of collected surfaces into different negative layers.")(
      "geo-tgeo-clayer-r-split", po::value<read_range>()->default_value({}),
      "R-tolerances (if > 0.) that triggers splitting "
      " of collected surfaces into different central layers.")(
      "geo-tgeo-clayer-z-split", po::value<read_range>()->default_value({}),
      "Z-tolerances (if > 0.) that triggers splitting "
      " of collected surfaces into different central layers.")(
      "geo-tgeo-player-r-split", po::value<read_range>()->default_value({}),
      "R-tolerances (if > 0.) that triggers splitting "
      " of collected surfaces into different positive layers.")(
      "geo-tgeo-sfbin-z-tolerance",
      po::value<std::vector<Interval>>()->default_value({}),
      "Tolerance interval in z [mm] for automated surface binning.")(
      "geo-tgeo-sfbin-r-tolerance",
      po::value<std::vector<Interval>>()->default_value({}),
      "Tolerance interval in r [mm] for automated surface binninng.")(
      "geo-tgeo-sfbin-phi-tolerance",
      po::value<std::vector<Interval>>()->default_value({}),
      "Tolerance interval in phi [rad] for automated surface binning.")(
      "geo-tgeo-nmodule-axes", po::value<read_strings>()->default_value({}),
      "Axes definition for negative sensitive objects, odered "
      "along the series.")(
      "geo-tgeo-player-z-split", po::value<read_range>()->default_value({}),
      "Z-tolerances (if > 0.) that triggers splitting "
      " of collected surfaces into different positive layers.")(
      "geo-tgeo-cmodule-axes", po::value<read_strings>()->default_value({}),
      "Axes definition for central sensitive objects, odered along the "
      "series.")("geo-tgeo-pmodule-axes",
                 po::value<read_strings>()->default_value({}),
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
    const variable_map_t& vm) {
  std::vector<Acts::TGeoLayerBuilder::Config> detLayerConfigs;

  // General: subdetector naming
  read_strings subDetectors =
      vm["geo-detector-volume"].template as<read_strings>();
  double unitScalor = vm["geo-tgeo-unit-scalor"].template as<double>();

  // The number of layers, can be set 1 with splitting detection
  read_series nLayers = vm["geo-tgeo-nlayers"].template as<read_series>();
  read_series cLayers = vm["geo-tgeo-clayers"].template as<read_series>();
  read_series pLayers = vm["geo-tgeo-players"].template as<read_series>();

  std::array<read_series, 3> layers = {nLayers, cLayers, pLayers};

  std::array<size_t, 3> seriesSize = {nLayers.size(), cLayers.size(),
                                      pLayers.size()};

  read_series ringLayout = vm["geo-tgeo-ringlayout"].template as<read_series>();

  read_range ringTolerance =
      vm["geo-tgeo-ringtolerance"].template as<read_range>();

  // The layer names to parse for in the TGeo
  read_strings nVolumeName =
      vm["geo-tgeo-nvolume-name"].template as<read_strings>();
  read_strings cVolumeName =
      vm["geo-tgeo-cvolume-name"].template as<read_strings>();
  read_strings pVolumeName =
      vm["geo-tgeo-pvolume-name"].template as<read_strings>();

  std::array<read_strings, 3> volumeName = {nVolumeName, cVolumeName,
                                            pVolumeName};

  read_strings nSensitiveNames =
      vm["geo-tgeo-nmodule-name"].template as<read_strings>();
  read_strings cSensitiveNames =
      vm["geo-tgeo-cmodule-name"].template as<read_strings>();
  read_strings pSensitiveNames =
      vm["geo-tgeo-pmodule-name"].template as<read_strings>();

  // The sensitive names to parse for in the TGeo
  std::array<read_strings, 3> sensitiveNames = {
      nSensitiveNames, cSensitiveNames, pSensitiveNames};

  read_strings nSensitiveAxes =
      vm["geo-tgeo-nmodule-axes"].template as<read_strings>();
  read_strings cSensitiveAxes =
      vm["geo-tgeo-cmodule-axes"].template as<read_strings>();
  read_strings pSensitiveAxes =
      vm["geo-tgeo-pmodule-axes"].template as<read_strings>();

  std::array<read_strings, 3> sensitiveAxes = {nSensitiveAxes, cSensitiveAxes,
                                               pSensitiveAxes};

  // The parse radii in r
  std::vector<Interval> nRrange =
      vm["geo-tgeo-nlayer-r-range"].template as<std::vector<Interval>>();
  std::vector<Interval> cRrange =
      vm["geo-tgeo-clayer-r-range"].template as<std::vector<Interval>>();
  std::vector<Interval> pRrange =
      vm["geo-tgeo-player-r-range"].template as<std::vector<Interval>>();
  std::array<std::vector<Interval>, 3> rRange = {nRrange, cRrange, pRrange};

  // The parse ranges in z
  std::vector<Interval> nZrange =
      vm["geo-tgeo-nlayer-z-range"].template as<std::vector<Interval>>();
  std::vector<Interval> cZrange =
      vm["geo-tgeo-clayer-z-range"].template as<std::vector<Interval>>();
  std::vector<Interval> pZrange =
      vm["geo-tgeo-player-z-range"].template as<std::vector<Interval>>();
  std::array<std::vector<Interval>, 3> zRange = {nZrange, cZrange, pZrange};

  // The split tolerances in r and z
  read_range nLayerSplitR =
      vm["geo-tgeo-nlayer-r-split"].template as<read_range>();
  read_range nLayerSplitZ =
      vm["geo-tgeo-nlayer-z-split"].template as<read_range>();
  read_range cLayerSplitR =
      vm["geo-tgeo-clayer-r-split"].template as<read_range>();
  read_range cLayerSplitZ =
      vm["geo-tgeo-clayer-z-split"].template as<read_range>();
  read_range pLayerSplitR =
      vm["geo-tgeo-player-r-split"].template as<read_range>();
  read_range pLayerSplitZ =
      vm["geo-tgeo-player-z-split"].template as<read_range>();

  std::array<read_range, 3> splitTolR = {nLayerSplitR, cLayerSplitR,
                                         pLayerSplitR};

  std::array<read_range, 3> splitTolZ = {nLayerSplitZ, cLayerSplitZ,
                                         pLayerSplitZ};

  // Automated binning configuration
  std::vector<Interval> binToleranceR =
      vm["geo-tgeo-sfbin-r-tolerance"].template as<std::vector<Interval>>();
  std::vector<Interval> binToleranceZ =
      vm["geo-tgeo-sfbin-z-tolerance"].template as<std::vector<Interval>>();
  std::vector<Interval> binTolerancePhi =
      vm["geo-tgeo-sfbin-phi-tolerance"].template as<std::vector<Interval>>();

  // The maximum series and total counter for access of nonsplit layers
  size_t max_series = *std::max_element(seriesSize.begin(), seriesSize.end());
  std::array<size_t, 3> ti = {0, 0, 0};

  // If a beam pipe is present, shift the sub detector names by one
  // Create a beam pipe if configured to do so
  int idetaddon = 0;
  auto beamPipeParameters =
      vm["geo-tgeo-bp-parameters"].template as<read_range>();
  if (beamPipeParameters.size() > 2) {
    ++idetaddon;
  }

  // Split the sensor names if there are mulitple ones
  auto splitAtOr =
      [](const std::string& sensorNames) -> std::vector<std::string> {
    std::vector<std::string> sensors;
    std::istringstream feed(sensorNames);
    std::string split;
    while (getline(feed, split, '|')) {
      sensors.push_back(split);
    }
    return sensors;
  };

  // Prepare the TGeoLayerBuilder::Configs
  for (size_t idet = 0; idet < max_series; ++idet) {
    // Each detector needs a layer builder
    Acts::TGeoLayerBuilder::Config layerBuilderConfig;

    // Configuration of AutoBinning
    if (binTolerancePhi.size() == binToleranceZ.size() and
        binToleranceR.size() == binToleranceZ.size()) {
      layerBuilderConfig.autoSurfaceBinning = true;

      auto tolR = binToleranceR[idet];
      std::vector<std::pair<double, double>> binTolerances{(int)Acts::binValues,
                                                           {0., 0.}};
      binTolerances[Acts::binR] = {tolR.lower.value_or(0.),
                                   tolR.upper.value_or(0.)};

      auto tolZ = binToleranceZ[idet];
      binTolerances[Acts::binZ] = {tolZ.lower.value_or(0.),
                                   tolZ.upper.value_or(0.)};
      auto tolPhi = binTolerancePhi[idet];
      binTolerances[Acts::binPhi] = {tolPhi.lower.value_or(0.),
                                     tolPhi.upper.value_or(0.)};

      layerBuilderConfig.surfaceBinMatcher =
          Acts::SurfaceBinningMatcher(binTolerances);
    }

    // Loop over the | n | c | p | configuration
    for (unsigned int ncp = 0; ncp < 3; ++ncp) {
      // number of layers of this configuration
      unsigned int nl = layers[ncp].size() > idet ? layers[ncp][idet] : 0;
      for (unsigned int in = 0; in < nl; ++in, ++ti[ncp]) {
        // Create the layer config object and fill it
        Acts::TGeoLayerBuilder::LayerConfig lConfig;
        lConfig.volumeName = volumeName[ncp][ti[ncp]];
        lConfig.sensorNames = splitAtOr(sensitiveNames[ncp][ti[ncp]]);
        lConfig.localAxes = sensitiveAxes[ncp][ti[ncp]];

        // Fill the parsing restrictions in r
        auto rR = rRange[ncp];
        if (rR.size() > ti[ncp]) {
          double rMin = rR[ti[ncp]].lower.value_or(0.);
          double rMax =
              rR[ti[ncp]].upper.value_or(std::numeric_limits<double>::max());
          lConfig.parseRanges.push_back({Acts::binR, {rMin, rMax}});
        }
        // Fill the layer splitting parameters in r
        if (splitTolR[ncp].size() > ti[ncp]) {
          double rsp = splitTolR[ncp][ti[ncp]];
          if (rsp > 0.) {
            lConfig.splitConfigs.push_back({Acts::binR, rsp});
          }
        }
        // Fill the parsing restrictions in z
        auto zR = zRange[ncp];
        if (zR.size() > ti[ncp]) {
          double zMin =
              zR[ti[ncp]].lower.value_or(-std::numeric_limits<double>::max());
          double zMax =
              zR[ti[ncp]].upper.value_or(std::numeric_limits<double>::max());
          lConfig.parseRanges.push_back({Acts::binZ, {zMin, zMax}});
        }
        // Fill the layer splitting parameters in z
        if (splitTolZ[ncp].size() > ti[ncp]) {
          double zsp = splitTolZ[ncp][ti[ncp]];
          if (zsp > 0.) {
            lConfig.splitConfigs.push_back({Acts::binZ, zsp});
          }
        }
        layerBuilderConfig.layerConfigurations[ncp].push_back(lConfig);
      }
    }

    // Set the scale and the layer creator
    layerBuilderConfig.configurationName = subDetectors[idet + idetaddon];
    layerBuilderConfig.unit = unitScalor;

    // Now add it to the configs
    detLayerConfigs.push_back(layerBuilderConfig);
  }
  return detLayerConfigs;
}

}  // namespace Options
}  // namespace ActsExamples

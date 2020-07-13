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
#include <sstream>
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
      "geo-tgeo-nodesearch-debug", po::value<bool>()->default_value(false),
      "Run ultra-verbose screen output to check node search.")(
      "geo-tgeo-unitScalor", po::value<double>()->default_value(10.),
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
      "geo-tgeo-nlayernames", po::value<read_strings>()->default_value({}),
      "Name identifier for negative layer objects, odered along the series.")(
      "geo-tgeo-clayernames", po::value<read_strings>()->default_value({}),
      "Name identifier for central layer objects, odered along the series.")(
      "geo-tgeo-playernames", po::value<read_strings>()->default_value({}),
      "Name identifier for positive layer objects, odered along the series.")(
      "geo-tgeo-nmodulenames", po::value<read_strings>()->default_value({}),
      "Name identifier for negative sensitive objects, odered along the "
      "series.")("geo-tgeo-cmodulenames",
                 po::value<read_strings>()->default_value({}),
                 "Name identifier for central sensitive objects, odered "
                 "along the series.")(
      "geo-tgeo-pmodulenames", po::value<read_strings>()->default_value({}),
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
      "geo-tgeo-nmoduleaxes", po::value<read_strings>()->default_value({}),
      "Axes definition for negative sensitive objects, odered "
      "along the series.")(
      "geo-tgeo-player-z-split", po::value<read_range>()->default_value({}),
      "Z-tolerances (if > 0.) that triggers splitting "
      " of collected surfaces into different positive layers.")(
      "geo-tgeo-cmoduleaxes", po::value<read_strings>()->default_value({}),
      "Axes definition for central sensitive objects, odered along the "
      "series.")("geo-tgeo-pmoduleaxes",
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
  read_strings subdetectors =
      vm["geo-detector-volume"].template as<read_strings>();
  double unitScalor = vm["geo-tgeo-unitScalor"].template as<double>();

  // The number of layers, can be set 1 with splitting detection
  read_series nlayers = vm["geo-tgeo-nlayers"].template as<read_series>();
  read_series clayers = vm["geo-tgeo-clayers"].template as<read_series>();
  read_series players = vm["geo-tgeo-players"].template as<read_series>();

  std::array<read_series, 3> layers = {nlayers, clayers, players};

  std::array<size_t, 3> series_size = {nlayers.size(), clayers.size(),
                                       players.size()};

  read_series ringlayout = vm["geo-tgeo-ringlayout"].template as<read_series>();

  read_range ringtolerance =
      vm["geo-tgeo-ringtolerance"].template as<read_range>();

  // The layer names to parse for in the TGeo
  read_strings nlayernames =
      vm["geo-tgeo-nlayernames"].template as<read_strings>();
  read_strings clayernames =
      vm["geo-tgeo-clayernames"].template as<read_strings>();
  read_strings playernames =
      vm["geo-tgeo-playernames"].template as<read_strings>();

  std::array<read_strings, 3> layernames = {nlayernames, clayernames,
                                            playernames};

  read_strings nsensitivenames =
      vm["geo-tgeo-nmodulenames"].template as<read_strings>();
  read_strings csensitivenames =
      vm["geo-tgeo-cmodulenames"].template as<read_strings>();
  read_strings psensitivenames =
      vm["geo-tgeo-pmodulenames"].template as<read_strings>();

  // The sensitive names to parse for in the TGeo
  std::array<read_strings, 3> sensitivenames = {
      nsensitivenames, csensitivenames, psensitivenames};

  read_strings nsensitiveaxes =
      vm["geo-tgeo-nmoduleaxes"].template as<read_strings>();
  read_strings csensitiveaxes =
      vm["geo-tgeo-cmoduleaxes"].template as<read_strings>();
  read_strings psensitiveaxes =
      vm["geo-tgeo-pmoduleaxes"].template as<read_strings>();

  std::array<read_strings, 3> sensitiveaxes = {nsensitiveaxes, csensitiveaxes,
                                               psensitiveaxes};

  // The parse radii in r
  std::vector<Interval> nrrange =
      vm["geo-tgeo-nlayer-r-range"].template as<std::vector<Interval>>();
  std::vector<Interval> crrange =
      vm["geo-tgeo-clayer-r-range"].template as<std::vector<Interval>>();
  std::vector<Interval> prrange =
      vm["geo-tgeo-player-r-range"].template as<std::vector<Interval>>();
  std::array<std::vector<Interval>, 3> rranges = {nrrange, crrange, prrange};

  // The parse ranges in z
  std::vector<Interval> nzrange =
      vm["geo-tgeo-nlayer-z-range"].template as<std::vector<Interval>>();
  std::vector<Interval> czrange =
      vm["geo-tgeo-clayer-z-range"].template as<std::vector<Interval>>();
  std::vector<Interval> pzrange =
      vm["geo-tgeo-player-z-range"].template as<std::vector<Interval>>();
  std::array<std::vector<Interval>, 3> zranges = {nzrange, czrange, pzrange};

  // The split tolerances in r and z
  read_range nlayersplitr =
      vm["geo-tgeo-nlayer-r-split"].template as<read_range>();
  read_range nlayersplitz =
      vm["geo-tgeo-nlayer-z-split"].template as<read_range>();
  read_range clayersplitr =
      vm["geo-tgeo-clayer-r-split"].template as<read_range>();
  read_range clayersplitz =
      vm["geo-tgeo-clayer-z-split"].template as<read_range>();
  read_range playersplitr =
      vm["geo-tgeo-player-r-split"].template as<read_range>();
  read_range playersplitz =
      vm["geo-tgeo-player-z-split"].template as<read_range>();

  std::array<read_range, 3> splittolr = {nlayersplitr, clayersplitr,
                                         playersplitr};

  std::array<read_range, 3> splittolz = {nlayersplitz, clayersplitz,
                                         playersplitz};

  // Automated binning configuration
  std::vector<Interval> binrtolerance =
      vm["geo-tgeo-sfbin-r-tolerance"].template as<std::vector<Interval>>();
  std::vector<Interval> binztolerance =
      vm["geo-tgeo-sfbin-z-tolerance"].template as<std::vector<Interval>>();
  std::vector<Interval> binphitolerance =
      vm["geo-tgeo-sfbin-phi-tolerance"].template as<std::vector<Interval>>();

  // The maximum series and total counter for access of nonsplit layers
  size_t max_series = *std::max_element(series_size.begin(), series_size.end());
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
    if (binphitolerance.size() == binrtolerance.size() and
        binrtolerance.size() == binztolerance.size()) {
      layerBuilderConfig.autoSurfaceBinning = true;

      auto rtol = binrtolerance[idet];
      std::vector<std::pair<double, double>> binTolerances{(int)Acts::binValues,
                                                           {0., 0.}};

      binTolerances[Acts::binR] = {rtol.lower.value_or(0.),
                                   rtol.upper.value_or(0.)};
      auto ztol = binztolerance[idet];
      binTolerances[Acts::binZ] = {ztol.lower.value_or(0.),
                                   ztol.upper.value_or(0.)};
      auto phitol = binphitolerance[idet];
      binTolerances[Acts::binPhi] = {phitol.lower.value_or(0.),
                                     phitol.upper.value_or(0.)};

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
        lConfig.layerName = layernames[ncp][ti[ncp]];
        lConfig.sensorNames = splitAtOr(sensitivenames[ncp][ti[ncp]]);
        lConfig.localAxes = sensitiveaxes[ncp][ti[ncp]];
        // Fill the parsing restrictions in r
        auto rrange = rranges[ncp];
        if (rrange.size() > ti[ncp]) {
          double rmin = rrange[ti[ncp]].lower.value_or(0.);
          double rmax = rrange[ti[ncp]].upper.value_or(
              std::numeric_limits<double>::max());
          lConfig.parseRanges.push_back({Acts::binR, {rmin, rmax}});
        }
        // Fill the layer splitting parameters in r
        if (splittolr[ncp].size() > ti[ncp]) {
          double rsplitol = splittolr[ncp][ti[ncp]];
          if (rsplitol > 0.) {
            lConfig.splitConfigs.push_back({Acts::binR, rsplitol});
          }
        }
        // Fill the parsing restrictions in z
        auto zrange = zranges[ncp];
        if (zrange.size() > ti[ncp]) {
          double zmin = zrange[ti[ncp]].lower.value_or(
              -std::numeric_limits<double>::max());
          double zmax = zrange[ti[ncp]].upper.value_or(
              std::numeric_limits<double>::max());
          lConfig.parseRanges.push_back({Acts::binZ, {zmin, zmax}});
        }
        // Fill the layer splitting parameters in z
        if (splittolz[ncp].size() > ti[ncp]) {
          double zsplitol = splittolz[ncp][ti[ncp]];
          if (zsplitol > 0.) {
            lConfig.splitConfigs.push_back({Acts::binZ, zsplitol});
          }
        }
        layerBuilderConfig.layerConfigurations[ncp].push_back(lConfig);
      }
    }

    // Set the scale and the layer creator
    layerBuilderConfig.configurationName = subdetectors[idet + idetaddon];
    layerBuilderConfig.unit = unitScalor;

    // Node search (ultra verbose) screen debug
    layerBuilderConfig.nodeSearchDebug =
        vm["geo-tgeo-nodesearch-debug"].template as<bool>();

    // Now add it to the configs
    detLayerConfigs.push_back(layerBuilderConfig);
  }
  return detLayerConfigs;
}

}  // namespace Options
}  // namespace FW

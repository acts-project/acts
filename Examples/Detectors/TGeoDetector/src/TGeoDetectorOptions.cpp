// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TGeoDetector/TGeoDetectorOptions.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <cstdlib>
#include <vector>

#include <boost/program_options.hpp>

void ActsExamples::Options::addTGeoGeometryOptions(Description& desc) {
  using boost::program_options::value;

  // due to the way program options handles options that can occur multiple
  // times, all options of a logical block must always be present; there are no
  // optional options.
  //
  // each detector volume configuration is one logical block which can
  // be repeated as many times as there are usable detector volumes.
  //
  // required per-volume options:
  //
  //   --geo-detector-volume InnerPixels
  //   --geo-tgeo-sfbin-r-tolerance 5:5
  //   --geo-tgeo-sfbin-phi-tolerance 0.025:0.025
  //   --geo-tgeo-sfbin-z-tolerance 5:5
  //
  // for each volume, a configurable number of negative/central/positive layers
  // can be configured. all `--geo-tgeo-{n,c,p}layers` option must always be
  // present and must be followed by the logical blocks for each layer. there
  // must be as many logical layer blocks as given in the corresponding
  // `*layers` option.
  //
  //   # two negative layers
  //   --geo-tgeo-nlayers 2
  //
  //   --geo-tgeo-nvolume-name Pixel::Pixel
  //   --geo-tgeo-nmodule-name Pixel::siLog
  //   --geo-tgeo-nmodule-axes YZX
  //   --geo-tgeo-nlayer-r-range 0:135
  //   --geo-tgeo-nlayer-z-range -3000:-250
  //   --geo-tgeo-nlayer-r-split -1.
  //   --geo-tgeo-nlayer-z-split 10.
  //
  //   --geo-tgeo-nvolume-name Pixel::SomethingElse
  //   --geo-tgeo-nmodule-name Pixel::FOO
  //   --geo-tgeo-nmodule-axes XYZ
  //   --geo-tgeo-nlayer-r-range 135:160
  //   --geo-tgeo-nlayer-z-range -3000:-250
  //   --geo-tgeo-nlayer-r-split -1.
  //   --geo-tgeo-nlayer-z-split 10.
  //
  //  # one central layer
  //   --geo-tgeo-clayers 1
  //
  //   --geo-tgeo-cvolume-name Pixel::Pixel
  //   --geo-tgeo-cmodule-name Pixel::siLog
  //   --geo-tgeo-cmodule-axes YZX
  //   --geo-tgeo-clayer-r-range 0:135
  //   --geo-tgeo-clayer-z-range -250:250
  //   --geo-tgeo-clayer-r-split 5.
  //   --geo-tgeo-clayer-z-split -1.
  //
  //   # no positive layer
  //   --geo-tgeo-players 0
  //
  //   # no --geo-tgeo-{cvolume,cmodule,clayer}* options
  //
  auto opt = desc.add_options();
  // required global options
  opt("geo-tgeo-filename", value<std::string>()->default_value(""),
      "Root file name.");
  opt("geo-tgeo-worldvolume", value<std::string>()->default_value(""),
      "Root world volume to start search from.")(
      "geo-tgeo-unit-scalor", value<double>()->default_value(10.),
      "Unit scalor from ROOT to Acts.");
  opt("geo-tgeo-bp-parameters", value<Doubles<3>>(),
      "Potential beam pipe parameters {r, z, t} in [mm].");
  // required per-volume options that can be present more than once
  opt("geo-tgeo-sfbin-r-tolerance", value<std::vector<Interval>>(),
      "Tolerance interval in r [mm] for automated surface binninng.");
  opt("geo-tgeo-sfbin-phi-tolerance", value<std::vector<Interval>>(),
      "Tolerance interval in phi [rad] for automated surface binning.");
  opt("geo-tgeo-sfbin-z-tolerance", value<std::vector<Interval>>(),
      "Tolerance interval in z [mm] for automated surface binning.");
  // optional per-volume layer options that can be present once.
  // `geo-tgeo-{n,c,p}-layers` must be present for each volume and if it is
  // non-zero, all other layer options with the same prefix must be present as
  // well.
  opt("geo-tgeo-nlayers", value<std::vector<int>>(),
      "Number of layers on the negative side.");
  opt("geo-tgeo-nvolume-name", value<std::vector<std::string>>(),
      "Name identifier of the volume for searching negative layers.");
  opt("geo-tgeo-nmodule-name", value<std::vector<std::string>>(),
      "Name identifier for negative sensitive objects.");
  opt("geo-tgeo-nmodule-axes", value<std::vector<std::string>>(),
      "Axes definition for negative sensitive objects.");
  opt("geo-tgeo-nlayer-r-range", value<std::vector<Interval>>(),
      "Radial range(s) for negative layers to restrict the module parsing.");
  opt("geo-tgeo-nlayer-z-range", value<std::vector<Interval>>(),
      "Longitudinal range(s) for negative layers to restrict the module "
      "parsing.");
  opt("geo-tgeo-nlayer-r-split", value<std::vector<double>>(),
      "R-tolerances (if > 0.) that triggers splitting "
      " of collected surfaces into different negative layers.");
  opt("geo-tgeo-nlayer-z-split", value<std::vector<double>>(),
      "Z-tolerances (if > 0.) that triggers splitting "
      " of collected surfaces into different negative layers.");
  // central layers options
  opt("geo-tgeo-clayers", value<std::vector<int>>(),
      "Number of layers in the barrel.");
  opt("geo-tgeo-cvolume-name", value<std::vector<std::string>>(),
      "Name identifier of the volume for searching central layers.");
  opt("geo-tgeo-cmodule-name", value<std::vector<std::string>>(),
      "Name identifier for central sensitive objects.");
  opt("geo-tgeo-cmodule-axes", value<std::vector<std::string>>(),
      "Axes definition for central sensitive objects.");
  opt("geo-tgeo-clayer-r-range", value<std::vector<Interval>>(),
      "Radial range(s) for central layers to restrict the module parsing.");
  opt("geo-tgeo-clayer-z-range", value<std::vector<Interval>>(),
      "Longitudinal range(s) for central layers to restrict the module "
      "parsing.");
  opt("geo-tgeo-clayer-r-split", value<std::vector<double>>(),
      "R-tolerances (if > 0.) that triggers splitting "
      " of collected surfaces into different central layers.");
  opt("geo-tgeo-clayer-z-split", value<std::vector<double>>(),
      "Z-tolerances (if > 0.) that triggers splitting "
      " of collected surfaces into different central layers.");
  // positive layers options
  opt("geo-tgeo-players", value<std::vector<int>>(),
      "Number of layers on the positive side.");
  opt("geo-tgeo-pvolume-name", value<std::vector<std::string>>(),
      "Name identifier of the volume for searching positive layers.");
  opt("geo-tgeo-pmodule-name", value<std::vector<std::string>>(),
      "Name identifier for positive sensitive objects.");
  opt("geo-tgeo-pmodule-axes", value<std::vector<std::string>>(),
      "Axes definition for positive sensitive objects.");
  opt("geo-tgeo-player-r-range", value<std::vector<Interval>>(),
      "Radial range(s) for positive layers to restrict the module parsing.");
  opt("geo-tgeo-player-z-range", value<std::vector<Interval>>(),
      "Longitudinal range(s) for positive layers to restrict the module "
      "parsing.");
  opt("geo-tgeo-player-r-split", value<std::vector<double>>(),
      "R-tolerances (if > 0.) that triggers splitting "
      " of collected surfaces into different positive layers.");
  opt("geo-tgeo-player-z-split", value<std::vector<double>>(),
      "Z-tolerances (if > 0.) that triggers splitting "
      " of collected surfaces into different positive layers.");
}

std::vector<Acts::TGeoLayerBuilder::Config>
ActsExamples::Options::readTGeoLayerBuilderConfigs(const Variables& vm) {
  std::vector<Acts::TGeoLayerBuilder::Config> detLayerConfigs;

  auto unitScalor = vm["geo-tgeo-unit-scalor"].template as<double>();

  // TODO how is this used?
  // If a beam pipe is present, shift the sub detector names by one
  //   size_t iVolOffset = vm.count("geo-tgeo-bp-parameters") ? 1 : 0;

  // subdetector selection
  auto subDetectors =
      vm["geo-detector-volume"].template as<std::vector<std::string>>();
  // per-volume automated binning configuration
  auto binToleranceR =
      vm["geo-tgeo-sfbin-r-tolerance"].template as<std::vector<Interval>>();
  auto binToleranceZ =
      vm["geo-tgeo-sfbin-z-tolerance"].template as<std::vector<Interval>>();
  auto binTolerancePhi =
      vm["geo-tgeo-sfbin-phi-tolerance"].template as<std::vector<Interval>>();
  // The number of layers, can be set -1 with automatic splitting detection
  std::array<std::vector<int>, 3> layers = {
      vm["geo-tgeo-nlayers"].template as<std::vector<int>>(),
      vm["geo-tgeo-clayers"].template as<std::vector<int>>(),
      vm["geo-tgeo-players"].template as<std::vector<int>>(),
  };
  // The layer names to parse for in the TGeo
  std::array<std::vector<std::string>, 3> volumeName = {
      vm["geo-tgeo-nvolume-name"].template as<std::vector<std::string>>(),
      vm["geo-tgeo-cvolume-name"].template as<std::vector<std::string>>(),
      vm["geo-tgeo-pvolume-name"].template as<std::vector<std::string>>(),
  };
  // The sensitive names to parse for in the TGeo
  std::array<std::vector<std::string>, 3> sensitiveNames = {
      vm["geo-tgeo-nmodule-name"].template as<std::vector<std::string>>(),
      vm["geo-tgeo-cmodule-name"].template as<std::vector<std::string>>(),
      vm["geo-tgeo-pmodule-name"].template as<std::vector<std::string>>(),
  };
  std::array<std::vector<std::string>, 3> sensitiveAxes = {
      vm["geo-tgeo-nmodule-axes"].template as<std::vector<std::string>>(),
      vm["geo-tgeo-cmodule-axes"].template as<std::vector<std::string>>(),
      vm["geo-tgeo-pmodule-axes"].template as<std::vector<std::string>>(),
  };
  // The parse radii in r
  std::array<std::vector<Interval>, 3> rRange = {
      vm["geo-tgeo-nlayer-r-range"].template as<std::vector<Interval>>(),
      vm["geo-tgeo-clayer-r-range"].template as<std::vector<Interval>>(),
      vm["geo-tgeo-player-r-range"].template as<std::vector<Interval>>()};
  // The parse ranges in z
  std::array<std::vector<Interval>, 3> zRange = {
      vm["geo-tgeo-nlayer-z-range"].template as<std::vector<Interval>>(),
      vm["geo-tgeo-clayer-z-range"].template as<std::vector<Interval>>(),
      vm["geo-tgeo-player-z-range"].template as<std::vector<Interval>>(),
  };
  // The split tolerances in r
  std::array<std::vector<double>, 3> splitTolR = {
      vm["geo-tgeo-nlayer-r-split"].template as<std::vector<double>>(),
      vm["geo-tgeo-clayer-r-split"].template as<std::vector<double>>(),
      vm["geo-tgeo-player-r-split"].template as<std::vector<double>>(),
  };
  // The split tolerances in r
  std::array<std::vector<double>, 3> splitTolZ = {
      vm["geo-tgeo-nlayer-z-split"].template as<std::vector<double>>(),
      vm["geo-tgeo-clayer-z-split"].template as<std::vector<double>>(),
      vm["geo-tgeo-player-z-split"].template as<std::vector<double>>(),
  };

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

  // iterate over all configured detector volumes
  //
  // current index within the negative/central/positive layers configuration
  // this has to be tracked separately as different volumes can have different
  // number of layers (or none at all) for each n/c/p configuration.
  std::array<size_t, 3> iLayers = {0, 0, 0};
  for (size_t iVol = 0; iVol < subDetectors.size(); ++iVol) {
    Acts::TGeoLayerBuilder::Config layerBuilderConfig;

    layerBuilderConfig.configurationName = subDetectors.at(iVol);
    layerBuilderConfig.unit = unitScalor;

    // configure surface autobinning
    // TODO how could this be optional
    auto tolR = binToleranceR.at(iVol);
    auto tolPhi = binTolerancePhi.at(iVol);
    auto tolZ = binToleranceZ.at(iVol);
    std::vector<std::pair<double, double>> binTolerances(
        static_cast<size_t>(Acts::binValues), {0., 0.});
    binTolerances[Acts::binR] = {tolR.lower.value_or(0.),
                                 tolR.upper.value_or(0.)};
    binTolerances[Acts::binZ] = {tolZ.lower.value_or(0.),
                                 tolZ.upper.value_or(0.)};
    binTolerances[Acts::binPhi] = {tolPhi.lower.value_or(0.),
                                   tolPhi.upper.value_or(0.)};
    layerBuilderConfig.autoSurfaceBinning = true;
    layerBuilderConfig.surfaceBinMatcher =
        Acts::SurfaceBinningMatcher(binTolerances);

    // loop over the negative/central/positive layer configurations
    for (size_t ncp = 0; ncp < 3; ++ncp) {
      // number of layers of this configuration
      size_t nLayers = layers[ncp].at(iVol);
      // location of the current volume layers in the configuration vectors
      size_t iLayer = iLayers[ncp];
      // WARNING do not use iLayers after this point in the iteration. it
      //   already contains the starting value for the next iteration.
      iLayers[ncp] += nLayers;

      // loop over all layers for this n/c/p volume configuration
      for (size_t endLayer = iLayer + nLayers; iLayer < endLayer; ++iLayer) {
        Acts::TGeoLayerBuilder::LayerConfig lConfig;

        lConfig.volumeName = volumeName[ncp].at(iLayer);
        lConfig.sensorNames = splitAtOr(sensitiveNames[ncp].at(iLayer));
        lConfig.localAxes = sensitiveAxes[ncp].at(iLayer);

        // Fill the parsing restrictions in r/z
        auto rR = rRange[ncp].at(iLayer);
        auto rMin = rR.lower.value_or(0.);
        auto rMax = rR.upper.value_or(std::numeric_limits<double>::max());
        auto zR = zRange[ncp].at(iLayer);
        auto zMin = zR.lower.value_or(-std::numeric_limits<double>::max());
        auto zMax = zR.upper.value_or(std::numeric_limits<double>::max());
        lConfig.parseRanges = {
            {Acts::binR, {rMin, rMax}},
            {Acts::binZ, {zMin, zMax}},
        };

        // Fill the layer splitting parameters in r/z
        auto str = splitTolR[ncp].at(iLayer);
        auto stz = splitTolZ[ncp].at(iLayer);
        if (0 < str) {
          lConfig.splitConfigs.emplace_back(Acts::binR, str);
        }
        if (0 < stz) {
          lConfig.splitConfigs.emplace_back(Acts::binZ, stz);
        }

        layerBuilderConfig.layerConfigurations[ncp].push_back(lConfig);
      }
    }

    // Now add it to the configs
    detLayerConfigs.push_back(layerBuilderConfig);
  }
  return detLayerConfigs;
}

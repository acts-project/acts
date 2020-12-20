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
  // times, all options of a logical block must always be present.
  //
  // each detector volume configuration is one logical block which can
  // be repeated as many times as there are usable detector volumes.
  //
  // required per-volume options:
  //
  //   --geo-tgeo-volume InnerPixels
  //   --geo-tgeo-sfbin-r-tolerance 5:5
  //   --geo-tgeo-sfbin-phi-tolerance 0.025:0.025
  //   --geo-tgeo-sfbin-z-tolerance 5:5
  //   --geo-tgeo-nlayers 0  # boolean switch whether there are negative layers
  //   --geo-tgeo-clayers 1  # boolean switch whether there are central layers
  //   --geo-tgeo-players 0  # boolean switch whether there are positive layers
  //
  // within each volume there can be negative/central/positive layers depending
  // on the which `--geo-tgeo-{n,c,p}layers` flags are set to true. if any of
  // them are set, they must be followed by the corresponding layer option. if
  // the `*layers` option is false, **no** further options **must** be set.
  //
  // examples: negative and central layers, but not positive layers
  //
  //   --geo-tgeo-nlayers 1
  //   --geo-tgeo-nvolume-name Pixel::Pixel
  //   --geo-tgeo-nmodule-name Pixel::siLog
  //   --geo-tgeo-nmodule-axes YZX
  //   --geo-tgeo-nlayer-r-range 0:135
  //   --geo-tgeo-nlayer-z-range -3000:-250
  //   --geo-tgeo-nlayer-r-split -1.
  //   --geo-tgeo-nlayer-z-split 10.
  //
  //   --geo-tgeo-clayers 1
  //   --geo-tgeo-cvolume-name Pixel::Pixel
  //   --geo-tgeo-cmodule-name Pixel::siLog
  //   --geo-tgeo-cmodule-axes YZX
  //   --geo-tgeo-clayer-r-range 0:135
  //   --geo-tgeo-clayer-z-range -250:250
  //   --geo-tgeo-clayer-r-split 5.
  //   --geo-tgeo-clayer-z-split -1.
  //
  //   --geo-tgeo-players 0
  //   # no --geo-tgeo-{cvolume,cmodule,clayer}* options
  //
  auto opt = desc.add_options();
  // required global options
  opt("geo-tgeo-filename", value<std::string>()->default_value(""),
      "Root file name.");
  opt("geo-tgeo-worldvolume", value<std::string>()->default_value(""),
      "Root world volume to start search from.");
  opt("geo-tgeo-unit-scalor", value<double>()->default_value(10.),
      "Unit scalor from ROOT to Acts.");
  opt("geo-tgeo-beampipe-parameters", value<Reals<3>>(),
      "Beam pipe parameters {r, z, t} in [mm]. Beam pipe is automatically "
      "created if the parameters are present.");
  // required per-volume options that can be present more than once
  opt("geo-tgeo-volume", value<std::vector<std::string>>(),
      "Detector volume name");
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
  opt("geo-tgeo-nlayers", value<std::vector<bool>>(),
      "Whether there are layers on the negative side.");
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
  opt("geo-tgeo-clayers", value<std::vector<bool>>(),
      "Whether there are layers in the central section.");
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
  opt("geo-tgeo-players", value<std::vector<bool>>(),
      "Whether there are layers on the positive side.");
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

  // subdetector selection
  auto subDetectors =
      vm["geo-tgeo-volume"].template as<std::vector<std::string>>();
  // per-volume automated binning configuration
  auto binToleranceR =
      vm["geo-tgeo-sfbin-r-tolerance"].template as<std::vector<Interval>>();
  auto binToleranceZ =
      vm["geo-tgeo-sfbin-z-tolerance"].template as<std::vector<Interval>>();
  auto binTolerancePhi =
      vm["geo-tgeo-sfbin-phi-tolerance"].template as<std::vector<Interval>>();
  // Whether any layers should be configured for a volume
  std::array<std::vector<bool>, 3> layers = {
      vm["geo-tgeo-nlayers"].template as<std::vector<bool>>(),
      vm["geo-tgeo-clayers"].template as<std::vector<bool>>(),
      vm["geo-tgeo-players"].template as<std::vector<bool>>(),
  };
  // The volume names to parse layers from in the TGeo
  std::array<std::vector<std::string>, 3> volumeName = {
      vm["geo-tgeo-nvolume-name"].template as<std::vector<std::string>>(),
      vm["geo-tgeo-cvolume-name"].template as<std::vector<std::string>>(),
      vm["geo-tgeo-pvolume-name"].template as<std::vector<std::string>>(),
  };
  // The sensitive surface/module names to parse for in the TGeo
  std::array<std::vector<std::string>, 3> sensitiveNames = {
      vm["geo-tgeo-nmodule-name"].template as<std::vector<std::string>>(),
      vm["geo-tgeo-cmodule-name"].template as<std::vector<std::string>>(),
      vm["geo-tgeo-pmodule-name"].template as<std::vector<std::string>>(),
  };
  // The sensitive surface axes configuration to parse for in the TGeo
  std::array<std::vector<std::string>, 3> sensitiveAxes = {
      vm["geo-tgeo-nmodule-axes"].template as<std::vector<std::string>>(),
      vm["geo-tgeo-cmodule-axes"].template as<std::vector<std::string>>(),
      vm["geo-tgeo-pmodule-axes"].template as<std::vector<std::string>>(),
  };
  // The layer transverse radius range
  std::array<std::vector<Interval>, 3> rRange = {
      vm["geo-tgeo-nlayer-r-range"].template as<std::vector<Interval>>(),
      vm["geo-tgeo-clayer-r-range"].template as<std::vector<Interval>>(),
      vm["geo-tgeo-player-r-range"].template as<std::vector<Interval>>()};
  // The layer z range
  std::array<std::vector<Interval>, 3> zRange = {
      vm["geo-tgeo-nlayer-z-range"].template as<std::vector<Interval>>(),
      vm["geo-tgeo-clayer-z-range"].template as<std::vector<Interval>>(),
      vm["geo-tgeo-player-z-range"].template as<std::vector<Interval>>(),
  };
  // The split tolerances in transverse radius
  std::array<std::vector<double>, 3> splitTolR = {
      vm["geo-tgeo-nlayer-r-split"].template as<std::vector<double>>(),
      vm["geo-tgeo-clayer-r-split"].template as<std::vector<double>>(),
      vm["geo-tgeo-player-r-split"].template as<std::vector<double>>(),
  };
  // The split tolerances in z
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
      if (not layers[ncp].at(iVol)) {
        continue;
      }

      // layer config position in the configuration vectors
      size_t iLayer = iLayers[ncp]++;

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

    detLayerConfigs.push_back(layerBuilderConfig);
  }
  return detLayerConfigs;
}

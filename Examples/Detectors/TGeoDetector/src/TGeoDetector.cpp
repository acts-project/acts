// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TGeoDetector/TGeoDetector.hpp"

#include "Acts/Geometry/CylinderVolumeBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/PassiveLayerBuilder.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Plugins/TGeo/TGeoCylinderDiscSplitter.hpp"
#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"
#include "Acts/Plugins/TGeo/TGeoLayerBuilder.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <list>

#include <boost/program_options.hpp>

namespace ActsExamples {
using namespace Options;

namespace {
/// @brief Function that constructs a set of layer builder
///        configs from a central @c TGeoDetector config.
///
/// @param config The input config
/// @return Vector of layer builder configs
std::vector<Acts::TGeoLayerBuilder::Config> makeLayerBuilderConfigs(
    const TGeoDetector::Config& config) {
  std::vector<Acts::TGeoLayerBuilder::Config> detLayerConfigs;

  // iterate over all configured detector volumes
  for (const auto& volume : config.volumes) {
    Acts::TGeoLayerBuilder::Config layerBuilderConfig;
    layerBuilderConfig.configurationName = volume.name;
    layerBuilderConfig.unit = config.unitScalor;

    // configure surface autobinning
    std::vector<std::pair<double, double>> binTolerances(
        static_cast<size_t>(Acts::binValues), {0., 0.});
    binTolerances[Acts::binR] = {volume.binToleranceR.lower.value_or(0.),
                                 volume.binToleranceR.upper.value_or(0.)};
    binTolerances[Acts::binZ] = {volume.binToleranceZ.lower.value_or(0.),
                                 volume.binToleranceZ.upper.value_or(0.)};
    binTolerances[Acts::binPhi] = {volume.binTolerancePhi.lower.value_or(0.),
                                   volume.binTolerancePhi.upper.value_or(0.)};
    layerBuilderConfig.autoSurfaceBinning = true;
    layerBuilderConfig.surfaceBinMatcher =
        Acts::SurfaceBinningMatcher(binTolerances);

    // loop over the negative/central/positive layer configurations
    for (auto ncp : {
             TGeoDetector::Config::Negative,
             TGeoDetector::Config::Central,
             TGeoDetector::Config::Positive,
         }) {
      if (!volume.layers.at(ncp)) {
        continue;
      }

      Acts::TGeoLayerBuilder::LayerConfig lConfig;
      lConfig.volumeName = volume.subVolumeName.at(ncp);
      lConfig.sensorNames = volume.sensitiveNames.at(ncp);
      lConfig.localAxes = volume.sensitiveAxes.at(ncp);

      // Fill the parsing restrictions in r/z
      auto rR = volume.rRange.at(ncp);
      auto rMin = rR.lower.value_or(0.);
      auto rMax = rR.upper.value_or(std::numeric_limits<double>::max());
      auto zR = volume.zRange.at(ncp);
      auto zMin = zR.lower.value_or(-std::numeric_limits<double>::max());
      auto zMax = zR.upper.value_or(std::numeric_limits<double>::max());
      lConfig.parseRanges = {
          {Acts::binR, {rMin, rMax}},
          {Acts::binZ, {zMin, zMax}},
      };

      // Fill the layer splitting parameters in r/z
      auto str = volume.splitTolR.at(ncp);
      auto stz = volume.splitTolZ.at(ncp);
      if (0 < str) {
        lConfig.splitConfigs.emplace_back(Acts::binR, str);
      }
      if (0 < stz) {
        lConfig.splitConfigs.emplace_back(Acts::binZ, stz);
      }

      layerBuilderConfig.layerConfigurations[ncp].push_back(lConfig);
    }

    // Perform splitting of cylinders and discs
    if (volume.cylinderDiscSplit) {
      Acts::TGeoCylinderDiscSplitter::Config cdsConfig;
      cdsConfig.cylinderPhiSegments = volume.cylinderNPhiSegments;
      cdsConfig.cylinderLongitudinalSegments = volume.cylinderNZSegments;
      cdsConfig.discPhiSegments = volume.discNPhiSegments;
      cdsConfig.discRadialSegments = volume.discNRSegments;
      layerBuilderConfig.detectorElementSplitter =
          std::make_shared<const Acts::TGeoCylinderDiscSplitter>(cdsConfig);
    }

    detLayerConfigs.push_back(layerBuilderConfig);
  }

  return detLayerConfigs;
}

/// @brief Function to build the generic tracking geometry
// from a TGeo object.
///
/// It does *currently* not translate the material, this has
/// to be done with a material mapping stage
///
/// @tparam variable_map_t is the variable map
///
/// @param vm is the variable map from the options
std::shared_ptr<const Acts::TrackingGeometry> buildTGeoDetector(
    const TGeoDetector::Config& config, const Acts::GeometryContext& context,
    std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>>&
        detElementStore,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator) {
  // configure surface array creator
  Acts::SurfaceArrayCreator::Config sacConfig;
  auto surfaceArrayCreator = std::make_shared<const Acts::SurfaceArrayCreator>(
      sacConfig,
      Acts::getDefaultLogger("SurfaceArrayCreator", config.surfaceLogLevel));
  // configure the proto layer helper
  Acts::ProtoLayerHelper::Config plhConfig;
  auto protoLayerHelper = std::make_shared<const Acts::ProtoLayerHelper>(
      plhConfig,
      Acts::getDefaultLogger("ProtoLayerHelper", config.layerLogLevel));
  // configure the layer creator that uses the surface array creator
  Acts::LayerCreator::Config lcConfig;
  lcConfig.surfaceArrayCreator = surfaceArrayCreator;
  auto layerCreator = std::make_shared<const Acts::LayerCreator>(
      lcConfig, Acts::getDefaultLogger("LayerCreator", config.layerLogLevel));
  // configure the layer array creator
  Acts::LayerArrayCreator::Config lacConfig;
  auto layerArrayCreator = std::make_shared<const Acts::LayerArrayCreator>(
      lacConfig,
      Acts::getDefaultLogger("LayerArrayCreator", config.layerLogLevel));
  // tracking volume array creator
  Acts::TrackingVolumeArrayCreator::Config tvacConfig;
  auto tVolumeArrayCreator =
      std::make_shared<const Acts::TrackingVolumeArrayCreator>(
          tvacConfig, Acts::getDefaultLogger("TrackingVolumeArrayCreator",
                                             config.volumeLogLevel));
  // configure the cylinder volume helper
  Acts::CylinderVolumeHelper::Config cvhConfig;
  cvhConfig.layerArrayCreator = layerArrayCreator;
  cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
  auto cylinderVolumeHelper =
      std::make_shared<const Acts::CylinderVolumeHelper>(
          cvhConfig, Acts::getDefaultLogger("CylinderVolumeHelper",
                                            config.volumeLogLevel));

  //-------------------------------------------------------------------------------------
  // list the volume builders
  std::list<std::shared_ptr<const Acts::ITrackingVolumeBuilder>> volumeBuilders;

  // Create a beam pipe if configured to do so
  if (config.buildBeamPipe) {
    /// configure the beam pipe layer builder
    Acts::PassiveLayerBuilder::Config bplConfig;
    bplConfig.layerIdentification = "BeamPipe";
    bplConfig.centralLayerRadii = {config.beamPipeRadius};
    bplConfig.centralLayerHalflengthZ = {config.beamPipeHalflengthZ};
    bplConfig.centralLayerThickness = {config.beamPipeLayerThickness};
    auto beamPipeBuilder = std::make_shared<const Acts::PassiveLayerBuilder>(
        bplConfig,
        Acts::getDefaultLogger("BeamPipeLayerBuilder", config.layerLogLevel));
    // create the volume for the beam pipe
    Acts::CylinderVolumeBuilder::Config bpvConfig;
    bpvConfig.trackingVolumeHelper = cylinderVolumeHelper;
    bpvConfig.volumeName = "BeamPipe";
    bpvConfig.layerBuilder = beamPipeBuilder;
    bpvConfig.layerEnvelopeR = {1. * Acts::UnitConstants::mm,
                                1. * Acts::UnitConstants::mm};
    bpvConfig.buildToRadiusZero = true;
    auto beamPipeVolumeBuilder =
        std::make_shared<const Acts::CylinderVolumeBuilder>(
            bpvConfig, Acts::getDefaultLogger("BeamPipeVolumeBuilder",
                                              config.volumeLogLevel));
    // add to the list of builders
    volumeBuilders.push_back(beamPipeVolumeBuilder);
  }

  // Import the file from
  TGeoManager::Import(config.fileName.c_str());

  auto layerBuilderConfigs = makeLayerBuilderConfigs(config);

  // Remember the layer builders to collect the detector elements
  std::vector<std::shared_ptr<const Acts::TGeoLayerBuilder>> tgLayerBuilders;

  for (auto& lbc : layerBuilderConfigs) {
    std::shared_ptr<const Acts::LayerCreator> layerCreatorLB = nullptr;

    if (lbc.autoSurfaceBinning) {
      // Configure surface array creator (optionally) per layer builder
      // (in order to configure them to work appropriately)
      Acts::SurfaceArrayCreator::Config sacConfigLB;
      sacConfigLB.surfaceMatcher = lbc.surfaceBinMatcher;
      auto surfaceArrayCreatorLB =
          std::make_shared<const Acts::SurfaceArrayCreator>(
              sacConfigLB, Acts::getDefaultLogger(
                               lbc.configurationName + "SurfaceArrayCreator",
                               config.surfaceLogLevel));
      // configure the layer creator that uses the surface array creator
      Acts::LayerCreator::Config lcConfigLB;
      lcConfigLB.surfaceArrayCreator = surfaceArrayCreatorLB;
      layerCreatorLB = std::make_shared<const Acts::LayerCreator>(
          lcConfigLB,
          Acts::getDefaultLogger(lbc.configurationName + "LayerCreator",
                                 config.layerLogLevel));
    }

    // Configure the proto layer helper
    Acts::ProtoLayerHelper::Config plhConfigLB;
    auto protoLayerHelperLB = std::make_shared<const Acts::ProtoLayerHelper>(
        plhConfigLB,
        Acts::getDefaultLogger(lbc.configurationName + "ProtoLayerHelper",
                               config.layerLogLevel));

    //-------------------------------------------------------------------------------------
    lbc.layerCreator =
        (layerCreatorLB != nullptr) ? layerCreatorLB : layerCreator;
    lbc.protoLayerHelper =
        (protoLayerHelperLB != nullptr) ? protoLayerHelperLB : protoLayerHelper;

    auto layerBuilder = std::make_shared<const Acts::TGeoLayerBuilder>(
        lbc, Acts::getDefaultLogger(lbc.configurationName + "LayerBuilder",
                                    config.layerLogLevel));
    // remember the layer builder
    tgLayerBuilders.push_back(layerBuilder);

    // build the pixel volume
    Acts::CylinderVolumeBuilder::Config volumeConfig;
    volumeConfig.trackingVolumeHelper = cylinderVolumeHelper;
    volumeConfig.volumeName = lbc.configurationName;
    volumeConfig.buildToRadiusZero = (volumeBuilders.size() == 0);
    volumeConfig.layerEnvelopeR = {1. * Acts::UnitConstants::mm,
                                   5. * Acts::UnitConstants::mm};
    auto ringLayoutConfiguration =
        [&](const std::vector<Acts::TGeoLayerBuilder::LayerConfig>& lConfigs)
        -> void {
      for (const auto& lcfg : lConfigs) {
        for (const auto& scfg : lcfg.splitConfigs) {
          if (scfg.first == Acts::binR and scfg.second > 0.) {
            volumeConfig.ringTolerance =
                std::max(volumeConfig.ringTolerance, scfg.second);
            volumeConfig.checkRingLayout = true;
          }
        }
      }
    };
    ringLayoutConfiguration(lbc.layerConfigurations[0]);
    ringLayoutConfiguration(lbc.layerConfigurations[2]);
    volumeConfig.layerBuilder = layerBuilder;
    volumeConfig.volumeSignature = 0;
    auto volumeBuilder = std::make_shared<const Acts::CylinderVolumeBuilder>(
        volumeConfig,
        Acts::getDefaultLogger(lbc.configurationName + "VolumeBuilder",
                               config.volumeLogLevel));
    // add to the list of builders
    volumeBuilders.push_back(volumeBuilder);
  }

  //-------------------------------------------------------------------------------------
  // create the tracking geometry
  Acts::TrackingGeometryBuilder::Config tgConfig;
  // Add the builders
  tgConfig.materialDecorator = mdecorator;

  for (auto& vb : volumeBuilders) {
    tgConfig.trackingVolumeBuilders.push_back(
        [=](const auto& gcontext, const auto& inner, const auto&) {
          return vb->trackingVolume(gcontext, inner);
        });
  }
  // Add the helper
  tgConfig.trackingVolumeHelper = cylinderVolumeHelper;
  auto cylinderGeometryBuilder =
      std::make_shared<const Acts::TrackingGeometryBuilder>(
          tgConfig, Acts::getDefaultLogger("TrackerGeometryBuilder",
                                           config.volumeLogLevel));
  // get the geometry
  auto trackingGeometry = cylinderGeometryBuilder->trackingGeometry(context);
  // collect the detector element store
  for (auto& lBuilder : tgLayerBuilders) {
    auto detElements = lBuilder->detectorElements();
    detElementStore.insert(detElementStore.begin(), detElements.begin(),
                           detElements.end());
  }

  /// return the tracking geometry
  return trackingGeometry;
}

/// Read the TGeo layer builder configurations from the user configuration.
void readTGeoLayerBuilderConfigs(const Variables& vm,
                                 TGeoDetector::Config& config) {
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

  // Check if layer builders are suggested to split cylinder / disk modules
  bool cylDiscSplit = vm["geo-tgeo-cyl-disc-split"].as<bool>();

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

  std::vector<int> cylZsegements = {};
  std::vector<int> cylPhiSegements = {};
  std::vector<int> discRsegements = {};
  std::vector<int> discPhiSegements = {};
  if (cylDiscSplit) {
    cylZsegements = vm["geo-tgeo-cyl-nz-segs"].template as<std::vector<int>>();
    cylPhiSegements =
        vm["geo-tgeo-cyl-nphi-segs"].template as<std::vector<int>>();
    discRsegements =
        vm["geo-tgeo-disc-nr-segs"].template as<std::vector<int>>();
    discPhiSegements =
        vm["geo-tgeo-disc-nphi-segs"].template as<std::vector<int>>();
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

  config.unitScalor = unitScalor;

  // iterate over all configured detector volumes
  //
  // current index within the negative/central/positive layers configuration
  // this has to be tracked separately as different volumes can have different
  // number of layers (or none at all) for each n/c/p configuration.
  std::array<size_t, 3> iLayers = {0, 0, 0};
  for (size_t iVol = 0; iVol < subDetectors.size(); ++iVol) {
    auto& volume = config.volumes.emplace_back();
    volume.name = subDetectors.at(iVol);

    volume.binToleranceR = binToleranceR.at(iVol);
    volume.binTolerancePhi = binTolerancePhi.at(iVol);
    volume.binToleranceZ = binToleranceZ.at(iVol);

    // loop over the negative/central/positive layer configurations
    for (size_t ncp = 0; ncp < 3; ++ncp) {
      if (not layers[ncp].at(iVol)) {
        continue;
      }

      // layer config position in the configuration vectors
      size_t iLayer = iLayers[ncp]++;

      TGeoDetector::Config::SubVolume eNcp{ncp};

      volume.layers.at(eNcp) = true;
      volume.subVolumeName.at(eNcp) = volumeName[ncp].at(iLayer);
      volume.sensitiveNames.at(eNcp) =
          splitAtOr(sensitiveNames[ncp].at(iLayer));
      volume.sensitiveAxes.at(eNcp) = sensitiveAxes[ncp].at(iLayer);

      // Fill the parsing restrictions in r/z
      volume.rRange.at(eNcp) = rRange[ncp].at(iLayer);
      volume.zRange.at(eNcp) = zRange[ncp].at(iLayer);

      // Fill the layer splitting parameters in r/z
      volume.splitTolR.at(eNcp) = splitTolR[ncp].at(iLayer);
      volume.splitTolZ.at(eNcp) = splitTolZ[ncp].at(iLayer);
    }

    // Perform splitting of cylinders and discs
    if (cylDiscSplit) {
      volume.cylinderDiscSplit = true;
      volume.cylinderNPhiSegments = cylPhiSegements.at(iVol);
      volume.cylinderNZSegments = cylZsegements.at(iVol);
      volume.discNPhiSegments = discPhiSegements.at(iVol);
      volume.discNRSegments = discRsegements.at(iVol);
    }
  }
}

}  // namespace

namespace ActsExamples {

void TGeoDetector::addOptions(
    boost::program_options::options_description& desc) const {
  using boost::program_options::value;

  // due to the way program options handles options that can occur multiple
  // times, all options of a logical block must always be present.
  //
  //   --geo-tgeo-cyl-disc-split # boolean switch to split cylinder/disc
  //   elements
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
  //  In case cylinder / disc splitting is on:
  //
  //   --geo-tgeo-cyl-nz-segs # number of z segments for cylinder splitting
  //   --geo-tgeo-cyl-nphi-segs # number of phi segments for cylinder splitting
  //   --geo-tgeo-disc-nr-segs # number of r segments for disc splitting
  //   --geo-tgeo-disc-nphi-segs # number of phi segments for disc splitting
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
  opt("geo-tgeo-cyl-disc-split", boost::program_options::bool_switch(),
      "Switch cylindrical / disc TGeo elements.");
  // required per-volume options that can be present more than once
  opt("geo-tgeo-volume", value<std::vector<std::string>>(),
      "Detector volume name");
  opt("geo-tgeo-sfbin-r-tolerance", value<std::vector<Interval>>(),
      "Tolerance interval in r [mm] for automated surface binninng.");
  opt("geo-tgeo-sfbin-phi-tolerance", value<std::vector<Interval>>(),
      "Tolerance interval in phi [rad] for automated surface binning.");
  opt("geo-tgeo-sfbin-z-tolerance", value<std::vector<Interval>>(),
      "Tolerance interval in z [mm] for automated surface binning.");
  // required per-volume IF 'geo-tgeo-cyl-disc-split' is set
  opt("geo-tgeo-cyl-nz-segs", value<std::vector<int>>(),
      "Conditional on geo-tgeo-cyl-disc-split: "
      "number of z segments for cylinder splitting.");
  opt("geo-tgeo-cyl-nphi-segs", value<std::vector<int>>(),
      "Conditional on geo-tgeo-cyl-disc-split: "
      "number of phi segments for cylinder splitting.");
  opt("geo-tgeo-disc-nr-segs", value<std::vector<int>>(),
      "Conditional on geo-tgeo-cyl-disc-split: "
      "number of r segments for disc splitting.");
  opt("geo-tgeo-disc-nphi-segs", value<std::vector<int>>(),
      "Conditional on geo-tgeo-cyl-disc-split: "
      "number of phi segments for disc splitting.");

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

auto TGeoDetector::finalize(
    const boost::program_options::variables_map& vm,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  Config config;
  config.surfaceLogLevel =
      Acts::Logging::Level(vm["geo-surface-loglevel"].template as<size_t>());
  config.layerLogLevel =
      Acts::Logging::Level(vm["geo-layer-loglevel"].template as<size_t>());
  config.volumeLogLevel =
      Acts::Logging::Level(vm["geo-volume-loglevel"].template as<size_t>());

  config.fileName = vm["geo-tgeo-filename"].template as<std::string>();

  config.buildBeamPipe = vm.count("geo-tgeo-beampipe-parameters") > 0;
  if (config.buildBeamPipe) {
    auto beamPipeParameters =
        vm["geo-tgeo-beampipe-parameters"].template as<Options::Reals<3>>();
    config.beamPipeRadius = beamPipeParameters[0];
    config.beamPipeHalflengthZ = beamPipeParameters[1];
    config.beamPipeLayerThickness = beamPipeParameters[2];
  }

  readTGeoLayerBuilderConfigs(vm, config);

  return finalize(config, std::move(mdecorator));
}

auto TGeoDetector::finalize(
    const Config& cfg,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  Acts::GeometryContext tGeoContext;
  TrackingGeometryPtr tgeoTrackingGeometry =
      buildTGeoDetector(cfg, tGeoContext, detectorStore, mdecorator);

  ContextDecorators tgeoContextDeocrators = {};
  // Return the pair of geometry and empty decorators
  return std::make_pair<TrackingGeometryPtr, ContextDecorators>(
      std::move(tgeoTrackingGeometry), std::move(tgeoContextDeocrators));
}

}  // namespace ActsExamples
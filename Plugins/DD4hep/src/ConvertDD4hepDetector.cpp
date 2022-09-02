// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/ConvertDD4hepDetector.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/PassiveLayerBuilder.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepMaterial.hpp"
#include "Acts/Plugins/DD4hep/DD4hepLayerBuilder.hpp"
#include "Acts/Plugins/DD4hep/DD4hepVolumeBuilder.hpp"
#include "Acts/Utilities/BinningData.hpp"

#include <list>
#include <regex>
#include <stdexcept>

#include "DD4hep/DetType.h"
#include "DDRec/DetectorData.h"
#include "TGeoManager.h"

namespace Acts {

namespace {
struct DebugVisitor {
  std::string operator()(int value) { return std::to_string(value); }

  std::string operator()(double value) { return std::to_string(value); }

  std::string operator()(std::string value) { return value; }
};
}  // namespace

std::unique_ptr<const TrackingGeometry> convertDD4hepDetector(
    dd4hep::DetElement worldDetElement, Logging::Level loggingLevel,
    BinningType bTypePhi, BinningType bTypeR, BinningType bTypeZ,
    double layerEnvelopeR, double layerEnvelopeZ, double defaultLayerThickness,
    const std::function<void(std::vector<dd4hep::DetElement>& detectors)>&
        sortSubDetectors,
    const Acts::GeometryContext& gctx,
    std::shared_ptr<const IMaterialDecorator> matDecorator,
    GeometryIdentifierHook geometryIdentifierHook) {
  // create local logger for conversion
  ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("DD4hepConversion", loggingLevel));
  ACTS_INFO("Translating DD4hep geometry into Acts geometry");
  // get the sub detectors of the world detector e.g. beampipe, pixel detector,
  // strip detector
  std::vector<dd4hep::DetElement> subDetectors;
  // go through the detector hierarchies
  collectSubDetectors_dd4hep(worldDetElement, subDetectors,
                             LoggerWrapper{*logger.m_logger});
  ACTS_VERBOSE("Collected " << subDetectors.size() << " sub detectors");
  // sort to build detector from bottom to top
  sortSubDetectors(subDetectors);
  // the volume builders of the subdetectors
  std::list<std::shared_ptr<const ITrackingVolumeBuilder>> volumeBuilders;
  // the beam pipe volume builder needs special treatment and needs to be added
  // in the end (beampipe exceeds length of all other subdetectors)
  std::shared_ptr<const CylinderVolumeBuilder> beamPipeVolumeBuilder;
  // loop over the sub detectors
  for (auto& subDetector : subDetectors) {
    ACTS_INFO("Translating DD4hep sub detector: " << subDetector.name());

    const dd4hep::rec::VariantParameters* params =
        subDetector.extension<dd4hep::rec::VariantParameters>(false);

    if (params != nullptr) {
      ACTS_VERBOSE("VariantParameters from DD4hep:");
      for (const auto& [k, v] : params->variantParameters) {
        ACTS_VERBOSE("- " << k << ": "
                          << boost::apply_visitor(DebugVisitor{}, v));
      }
    }

    // create volume builder
    auto volBuilder = volumeBuilder_dd4hep(
        subDetector, loggingLevel, bTypePhi, bTypeR, bTypeZ, layerEnvelopeR,
        layerEnvelopeZ, defaultLayerThickness);
    if (volBuilder) {
      // distinguish beam pipe
      if (volBuilder->getConfiguration().buildToRadiusZero) {
        // check if beam pipe is already present
        if (beamPipeVolumeBuilder) {
          throw std::logic_error(
              std::string("Beampipe has already been set! There can only "
                          "exist one beam pipe. Please check your "
                          "detector construction. Current volume name: ") +
              volBuilder->getConfiguration().volumeName +
              std::string(", name of volume, already set as beam pipe: ") +
              beamPipeVolumeBuilder->getConfiguration().volumeName);
        }
        // set the beam pipe
        beamPipeVolumeBuilder = volBuilder;
      } else {
        volumeBuilders.push_back(volBuilder);
      }
    }
  }
  // Finally add the beam pipe
  if (beamPipeVolumeBuilder) {
    volumeBuilders.push_back(beamPipeVolumeBuilder);
  }

  std::vector<std::function<std::shared_ptr<TrackingVolume>(
      const GeometryContext&, const TrackingVolumePtr&,
      const VolumeBoundsPtr&)>>
      volumeFactories;

  for (const auto& vb : volumeBuilders) {
    volumeFactories.push_back(
        [vb](const GeometryContext& vgctx,
             const std::shared_ptr<const TrackingVolume>& inner,
             const VolumeBoundsPtr&) {
          return vb->trackingVolume(vgctx, inner);
        });
  }

  // create cylinder volume helper
  auto volumeHelper = cylinderVolumeHelper_dd4hep();
  // hand over the collected volume builders
  Acts::TrackingGeometryBuilder::Config tgbConfig;
  tgbConfig.trackingVolumeHelper = volumeHelper;
  tgbConfig.materialDecorator = std::move(matDecorator);
  tgbConfig.trackingVolumeBuilders = std::move(volumeFactories);
  tgbConfig.geometryIdentifierHook = std::move(geometryIdentifierHook);
  auto trackingGeometryBuilder =
      std::make_shared<const Acts::TrackingGeometryBuilder>(tgbConfig);
  return (trackingGeometryBuilder->trackingGeometry(gctx));
}

std::shared_ptr<const CylinderVolumeBuilder> volumeBuilder_dd4hep(
    dd4hep::DetElement subDetector, Logging::Level loggingLevel,
    BinningType bTypePhi, BinningType bTypeR, BinningType bTypeZ,
    double layerEnvelopeR, double layerEnvelopeZ,
    double defaultLayerThickness) {
  // create cylinder volume helper
  auto volumeHelper = cylinderVolumeHelper_dd4hep(loggingLevel);
  // create local logger for conversion
  ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("D2A_Logger", loggingLevel));
  ACTS_VERBOSE("Processing detector element:  " << subDetector.name());

  LoggerWrapper loggerWrapper{*logger.m_logger};

  dd4hep::DetType subDetType{subDetector.typeFlag()};
  ACTS_VERBOSE("SubDetector type is: ["
               << subDetType << "], compound: "
               << (subDetector.type() == "compound" ? "yes" : "no"));
  if (subDetector.type() == "compound") {
    ACTS_VERBOSE("Subdetector : '" << subDetector.name()
                                   << "' has type compound ");
    ACTS_VERBOSE(
        "handling as a compound volume (a hierachy of a "
        "barrel-endcap structure) and resolving the "
        "subvolumes...");
    // Now create the Layerbuilders and Volumebuilder
    // the layers
    /// the dd4hep::DetElements of the layers of the negative volume
    std::vector<dd4hep::DetElement> negativeLayers;
    /// the dd4hep::DetElements of the layers of the central volume
    std::vector<dd4hep::DetElement> centralLayers;
    /// the dd4hep::DetElements of the layers of the positive volume
    std::vector<dd4hep::DetElement> positiveLayers;

    // the configuration object of the volume builder
    Acts::CylinderVolumeBuilder::Config cvbConfig;

    // go through sub volumes
    std::vector<dd4hep::DetElement> compounds;
    collectCompounds_dd4hep(subDetector, compounds);

    // get z position to distinguish positive & negative endcap
    double zPos = 0.;
    // flags to catch if sub volumes have been set already
    bool nEndCap = false;
    bool pEndCap = false;
    bool barrel = false;
    for (auto& volumeDetElement : compounds) {
      ACTS_VERBOSE("Volume : '"
                   << subDetector.name()
                   << "' is a compound volume -> resolve the sub volumes");

      // get the dimensions of the volume
      TGeoShape* geoShape =
          volumeDetElement.placement().ptr()->GetVolume()->GetShape();
      // check if it has a shape (the other case should not happen)
      if (geoShape != nullptr) {
        zPos = volumeDetElement.placement()
                   .ptr()
                   ->GetMatrix()
                   ->GetTranslation()[2] *
               UnitConstants::cm;
      } else {
        throw std::logic_error(std::string("Volume of DetElement: ") +
                               volumeDetElement.name() +
                               std::string(" has no shape!"));
      }

      dd4hep::DetType type{volumeDetElement.typeFlag()};

      if (type.is(dd4hep::DetType::ENDCAP)) {
        ACTS_VERBOSE(std::string("Subvolume : '") + volumeDetElement.name() +
                     std::string("' is marked ENDCAP"));
        if (zPos < 0.) {
          if (nEndCap) {
            throw std::logic_error(
                "Negative Endcap was already given for this "
                "hierachy! Please create a new "
                "DD4hep_SubDetectorAssembly for the next "
                "hierarchy.");
          }
          nEndCap = true;
          ACTS_VERBOSE("-> is negative endcap");
          ACTS_VERBOSE("-> collecting layers");
          collectLayers_dd4hep(volumeDetElement, negativeLayers, loggerWrapper);
          // Fill the volume material for barrel case
          if (getParamOr<bool>("boundary_material", volumeDetElement, false)) {
            ACTS_VERBOSE(
                "-> boundary_material flag detected, creating proto "
                "material.");
            auto& params = getParams(volumeDetElement);
            if (hasParam("boundary_material_negative", volumeDetElement)) {
              ACTS_VERBOSE("--> negative");
              cvbConfig.boundaryMaterial[2] = Acts::createProtoMaterial(
                  params, "boundary_material_negative",
                  {{"binPhi", Acts::closed}, {"binR", Acts::open}},
                  LoggerWrapper{logger()});
            }
            if (hasParam("boundary_material_positive", volumeDetElement)) {
              ACTS_VERBOSE("--> positive");
              cvbConfig.boundaryMaterial[3] = Acts::createProtoMaterial(
                  params, "boundary_material_positive",
                  {{"binPhi", Acts::closed}, {"binR", Acts::open}},
                  LoggerWrapper{logger()});
            }
          }
        } else {
          if (pEndCap) {
            throw std::logic_error(
                "Positive Endcap was already given for this "
                "hierachy! Please create a new "
                "DD4hep_SubDetectorAssembly for the next "
                "hierarchy.");
          }
          pEndCap = true;
          ACTS_VERBOSE("-> is positive endcap");
          ACTS_VERBOSE("-> collecting layers");
          collectLayers_dd4hep(volumeDetElement, positiveLayers, loggerWrapper);
          // Fill the volume material for barrel case
          if (getParamOr<bool>("boundary_material", volumeDetElement, false)) {
            ACTS_VERBOSE(
                "-> boundary_material flag detected, creating proto "
                "material.");
            auto& params = getParams(volumeDetElement);
            if (params.contains("boundary_material_negative")) {
              ACTS_VERBOSE("--> negative");
              cvbConfig.boundaryMaterial[4] = Acts::createProtoMaterial(
                  params, "boundary_material_negative",
                  {{"binPhi", Acts::closed}, {"binR", Acts::open}},
                  LoggerWrapper{logger()});
            }
            if (params.contains("boundary_material_positive")) {
              ACTS_VERBOSE("--> positive");
              cvbConfig.boundaryMaterial[5] = Acts::createProtoMaterial(
                  params, "boundary_material_positive",
                  {{"binPhi", Acts::closed}, {"binR", Acts::open}},
                  LoggerWrapper{logger()});
            }
          }
        }
      } else if (type.is(dd4hep::DetType::BARREL)) {
        if (barrel) {
          throw std::logic_error(
              "Barrel was already given for this "
              "hierachy! Please create a new "
              "DD4hep_SubDetectorAssembly for the next "
              "hierarchy.");
        }
        barrel = true;
        ACTS_VERBOSE("Subvolume : " << volumeDetElement.name()
                                    << " is marked as BARREL");
        ACTS_VERBOSE("-> collecting layers");
        collectLayers_dd4hep(volumeDetElement, centralLayers, loggerWrapper);
        // Fill the volume material for barrel case
        if (getParamOr<bool>("boundary_material", volumeDetElement, false)) {
          ACTS_VERBOSE(
              "-> boundary_material flag detected, creating proto "
              "material.");
          auto& params = getParams(volumeDetElement);
          if (params.contains("boundary_material_negative")) {
            ACTS_VERBOSE("--> negative");
            cvbConfig.boundaryMaterial[3] = Acts::createProtoMaterial(
                params, "boundary_material_negative",
                {{"binPhi", Acts::closed}, {"binR", Acts::open}},
                LoggerWrapper{logger()});
          }
          if (params.contains("boundary_material_positive")) {
            ACTS_VERBOSE("--> positive");
            cvbConfig.boundaryMaterial[4] = Acts::createProtoMaterial(
                params, "boundary_material_positive",
                {{"binPhi", Acts::closed}, {"binR", Acts::open}},
                LoggerWrapper{logger()});
          }
        }
      } else {
        throw std::logic_error(
            std::string("Current DetElement: ") + volumeDetElement.name() +
            std::string(" has inconsistent settings. It's a compound,"
                        " but its DetectorType is neither BARREL nor ENDCAP"
                        " Please check your detector construction."));
      }

      // Fill the volume material for the inner / outer cover
      if (getParamOr<bool>("boundary_material", volumeDetElement, false)) {
        ACTS_VERBOSE(
            "-> boundary_material flag detected, creating proto "
            "material.");
        auto& params = getParams(volumeDetElement);
        if (params.contains("boundary_material_inner")) {
          ACTS_VERBOSE("--> inner");
          cvbConfig.boundaryMaterial[0] = Acts::createProtoMaterial(
              params, "boundary_material_inner",
              {{"binPhi", Acts::closed}, {"binZ", Acts::open}},
              LoggerWrapper{logger()});
        }
        if (params.contains("boundary_material_outer")) {
          ACTS_VERBOSE("--> outer");
          cvbConfig.boundaryMaterial[1] = Acts::createProtoMaterial(
              params, "boundary_material_outer",
              {{"binPhi", Acts::closed}, {"binZ", Acts::open}},
              LoggerWrapper{logger()});
        }
      }
    }

    if ((pEndCap && !nEndCap) || (!pEndCap && nEndCap)) {
      throw std::logic_error(
          "Only one Endcap is given for the current "
          "hierarchy! Endcaps should always occur in "
          "pairs. Please check your detector "
          "construction.");
    }

    // configure SurfaceArrayCreator
    auto surfaceArrayCreator =
        std::make_shared<const Acts::SurfaceArrayCreator>(
            Acts::getDefaultLogger("D2A_SAC", loggingLevel));
    // configure LayerCreator
    Acts::LayerCreator::Config lcConfig;
    lcConfig.surfaceArrayCreator = surfaceArrayCreator;
    auto layerCreator = std::make_shared<const Acts::LayerCreator>(
        lcConfig, Acts::getDefaultLogger("D2A_LAC", loggingLevel));
    // configure DD4hepLayerBuilder
    Acts::DD4hepLayerBuilder::Config lbConfig;
    lbConfig.configurationName = subDetector.name();
    lbConfig.layerCreator = layerCreator;
    lbConfig.negativeLayers = negativeLayers;
    lbConfig.centralLayers = centralLayers;
    lbConfig.positiveLayers = positiveLayers;
    lbConfig.bTypePhi = bTypePhi;
    lbConfig.bTypeR = bTypeR;
    lbConfig.bTypeZ = bTypeZ;
    lbConfig.defaultThickness = defaultLayerThickness;
    auto dd4hepLayerBuilder = std::make_shared<const Acts::DD4hepLayerBuilder>(
        lbConfig,
        Acts::getDefaultLogger(std::string("D2A_L:") + subDetector.name(),
                               loggingLevel));

    // Create the sub volume
    // Dimensions are created automatically by adding a tolerance to the
    // layer setup
    cvbConfig.layerEnvelopeR = std::make_pair(layerEnvelopeR, layerEnvelopeR);
    cvbConfig.layerEnvelopeZ = layerEnvelopeZ;
    cvbConfig.trackingVolumeHelper = volumeHelper;
    cvbConfig.volumeSignature = 0;
    cvbConfig.volumeName = subDetector.name();
    cvbConfig.layerBuilder = dd4hepLayerBuilder;
    auto cylinderVolumeBuilder =
        std::make_shared<const Acts::CylinderVolumeBuilder>(
            cvbConfig,
            Acts::getDefaultLogger(std::string("D2A_V:") + subDetector.name(),
                                   loggingLevel));
    return cylinderVolumeBuilder;
  } else if (subDetType.is(dd4hep::DetType::BEAMPIPE) ||
             getParamOr<bool>("passive_layer", subDetector, false)) {
    ACTS_VERBOSE("Subdetector : " << subDetector.name()
                                  << " - building a passive cylinder.");

    if (subDetType.is(dd4hep::DetType::BEAMPIPE)) {
      ACTS_VERBOSE("This is the beam pipe - will be built to r -> 0.");
    }

    // get the dimensions of the volume
    TGeoShape* geoShape =
        subDetector.placement().ptr()->GetVolume()->GetShape();
    TGeoTubeSeg* tube = dynamic_cast<TGeoTubeSeg*>(geoShape);
    if (tube == nullptr) {
      throw std::logic_error(
          "Cylinder has wrong shape - needs to be TGeoTubeSeg!");
    }
    // get the dimension of TGeo and convert lengths
    double rMin = tube->GetRmin() * UnitConstants::cm - layerEnvelopeR;
    double rMax = tube->GetRmax() * UnitConstants::cm + layerEnvelopeR;
    double halfZ = tube->GetDz() * UnitConstants::cm + layerEnvelopeZ;
    ACTS_VERBOSE(
        "Extracting cylindrical volume bounds ( rmin / rmax / "
        "halfZ )=  ( "
        << rMin << " / " << rMax << " / " << halfZ << " )");

    std::shared_ptr<Acts::ISurfaceMaterial> plMaterial = nullptr;
    if (getParamOr<bool>("layer_material", subDetector, false)) {
      // get the possible material of the surounding volume
      ACTS_VERBOSE("--> adding layer material at 'representing'");
      plMaterial = Acts::createProtoMaterial(
          getParams(subDetector), "layer_material_representing",
          {{"binPhi", Acts::closed}, {"binZ", Acts::open}},
          LoggerWrapper{logger()});
    }

    // configure the passive layer builder
    Acts::PassiveLayerBuilder::Config plbConfig;
    plbConfig.layerIdentification = subDetector.name();
    plbConfig.centralLayerRadii = std::vector<double>(1, 0.5 * (rMax + rMin));
    plbConfig.centralLayerHalflengthZ = std::vector<double>(1, halfZ);
    plbConfig.centralLayerThickness = std::vector<double>(1, fabs(rMax - rMin));
    plbConfig.centralLayerMaterial = {plMaterial};
    auto pcLayerBuilder = std::make_shared<const Acts::PassiveLayerBuilder>(
        plbConfig,
        Acts::getDefaultLogger(std::string("D2A_PL:") + subDetector.name(),
                               loggingLevel));

    // the configuration object of the volume builder
    Acts::CylinderVolumeBuilder::Config cvbConfig;
    cvbConfig.trackingVolumeHelper = volumeHelper;
    cvbConfig.volumeSignature = 0;
    cvbConfig.volumeName = subDetector.name();
    cvbConfig.layerBuilder = pcLayerBuilder;
    cvbConfig.layerEnvelopeR = {layerEnvelopeR, layerEnvelopeR};
    cvbConfig.layerEnvelopeZ = layerEnvelopeZ;
    cvbConfig.buildToRadiusZero = subDetType.is(dd4hep::DetType::BEAMPIPE);

    // Fill the volume material for the inner / outer cover
    if (getParamOr<bool>("boundary_material", subDetector, false)) {
      ACTS_VERBOSE(
          "-> boundary_material flag detected, creating proto "
          "material.");
      auto& params = getParams(subDetector);
      if (hasParam("boundary_material_inner", subDetector)) {
        ACTS_VERBOSE("--> inner");
        cvbConfig.boundaryMaterial[0] = Acts::createProtoMaterial(
            params, "boundary_material_inner",
            {{"binPhi", Acts::closed}, {"binZ", Acts::open}},
            LoggerWrapper{logger()});
      }
      if (hasParam("boundary_material_outer", subDetector)) {
        ACTS_VERBOSE("--> outer");
        cvbConfig.boundaryMaterial[1] = Acts::createProtoMaterial(
            params, "boundary_material_outer",
            {{"binPhi", Acts::closed}, {"binZ", Acts::open}},
            LoggerWrapper{logger()});
      }
    }

    // beam pipe / passive cylinder volume builder
    auto pcVolumeBuilder = std::make_shared<const Acts::CylinderVolumeBuilder>(
        cvbConfig,
        Acts::getDefaultLogger(std::string("D2A_V:") + subDetector.name(),
                               loggingLevel));
    return pcVolumeBuilder;
  } else if (subDetType.is(dd4hep::DetType::BARREL)) {
    ACTS_VERBOSE("Subdetector: "
                 << subDetector.name()
                 << " is a (sensitive) Barrel volume - building barrel.");
    /// the dd4hep::DetElements of the layers of the central volume
    std::vector<dd4hep::DetElement> centralLayers, centralVolumes;
    ACTS_VERBOSE("-> collecting layers");
    collectLayers_dd4hep(subDetector, centralLayers, loggerWrapper);

    // configure SurfaceArrayCreator
    auto surfaceArrayCreator =
        std::make_shared<const Acts::SurfaceArrayCreator>(
            Acts::getDefaultLogger("D2A_SAC", loggingLevel));
    // configure LayerCreator
    Acts::LayerCreator::Config lcConfig;
    lcConfig.surfaceArrayCreator = surfaceArrayCreator;
    auto layerCreator = std::make_shared<const Acts::LayerCreator>(
        lcConfig, Acts::getDefaultLogger("D2A_LAC", loggingLevel));
    // configure DD4hepLayerBuilder
    Acts::DD4hepLayerBuilder::Config lbConfig;
    lbConfig.configurationName = subDetector.name();
    lbConfig.layerCreator = layerCreator;
    lbConfig.centralLayers = centralLayers;
    lbConfig.bTypePhi = bTypePhi;
    lbConfig.bTypeZ = bTypeZ;
    lbConfig.defaultThickness = defaultLayerThickness;
    auto dd4hepLayerBuilder = std::make_shared<const Acts::DD4hepLayerBuilder>(
        lbConfig,
        Acts::getDefaultLogger(std::string("D2A_LB_") + subDetector.name(),
                               loggingLevel));

    // Configure DD4hepVolumeBuilder
    Acts::DD4hepVolumeBuilder::Config vbConfig;
    vbConfig.configurationName = subDetector.name();
    vbConfig.centralVolumes = centralVolumes;
    auto dd4hepVolumeBuilder =
        std::make_shared<const Acts::DD4hepVolumeBuilder>(
            vbConfig,
            Acts::getDefaultLogger(std::string("D2A_VB_") + subDetector.name(),
                                   loggingLevel));

    // the configuration object of the volume builder
    Acts::CylinderVolumeBuilder::Config cvbConfig;
    // get the dimensions of the volume
    TGeoShape* geoShape =
        subDetector.placement().ptr()->GetVolume()->GetShape();
    // this should not happen
    if (geoShape == nullptr) {
      throw std::logic_error(std::string("Volume of DetElement: ") +
                             subDetector.name() +
                             std::string(" has no a shape!"));
    }

    cvbConfig.layerEnvelopeR = std::make_pair(layerEnvelopeR, layerEnvelopeR);
    cvbConfig.layerEnvelopeZ = layerEnvelopeZ;
    cvbConfig.trackingVolumeHelper = volumeHelper;
    cvbConfig.volumeSignature = 0;
    cvbConfig.volumeName = subDetector.name();
    cvbConfig.layerBuilder = dd4hepLayerBuilder;
    cvbConfig.ctVolumeBuilder = dd4hepVolumeBuilder;
    auto cylinderVolumeBuilder =
        std::make_shared<const Acts::CylinderVolumeBuilder>(
            cvbConfig,
            Acts::getDefaultLogger(std::string("D2A_V:") + subDetector.name(),
                                   loggingLevel));
    return cylinderVolumeBuilder;
  } else {
    ACTS_INFO(
        "Subdetector with name : '"
        << subDetector.name()
        << "' has inconsistent information for translation and is not of type "
           "'compound'. If you want to have this DetElement be translated "
           "into the tracking geometry you need add the right DetectorType "
           "or VariantParameters (at this stage the subvolume needs to be "
           "declared as BEAMPIPE or BARREl, or have a VariantParameter "
           "passive_layer=true) or if it is a compound DetElement (containing "
           "a barrel-endcap hierarchy), the type needs to be set to "
           "'compound'.");
    return nullptr;
  }
}

std::shared_ptr<const Acts::CylinderVolumeHelper> cylinderVolumeHelper_dd4hep(
    Logging::Level loggingLevel) {
  // create cylindervolumehelper which can be used by all instances
  // hand over LayerArrayCreator
  Acts::LayerArrayCreator::Config lacConfig;
  auto layerArrayCreator = std::make_shared<const Acts::LayerArrayCreator>(
      lacConfig, Acts::getDefaultLogger(std::string("D2A_LAC"), loggingLevel));
  // tracking volume array creator
  Acts::TrackingVolumeArrayCreator::Config tvacConfig;
  auto trackingVolumeArrayCreator =
      std::make_shared<const Acts::TrackingVolumeArrayCreator>(
          tvacConfig,
          Acts::getDefaultLogger(std::string("D2A_TVAC"), loggingLevel));
  // configure the cylinder volume helper
  Acts::CylinderVolumeHelper::Config cvhConfig;
  cvhConfig.layerArrayCreator = layerArrayCreator;
  cvhConfig.trackingVolumeArrayCreator = trackingVolumeArrayCreator;
  auto cylinderVolumeHelper =
      std::make_shared<const Acts::CylinderVolumeHelper>(
          cvhConfig,
          Acts::getDefaultLogger(std::string("D2A_CVH"), loggingLevel));

  return cylinderVolumeHelper;
}

void collectCompounds_dd4hep(dd4hep::DetElement& detElement,
                             std::vector<dd4hep::DetElement>& compounds) {
  const dd4hep::DetElement::Children& children = detElement.children();
  for (auto& child : children) {
    dd4hep::DetElement childDetElement = child.second;
    dd4hep::DetType type{childDetElement.typeFlag()};
    if (type.is(dd4hep::DetType::BARREL) || type.is(dd4hep::DetType::ENDCAP)) {
      compounds.push_back(childDetElement);
    }
    collectCompounds_dd4hep(childDetElement, compounds);
  }
}

void collectSubDetectors_dd4hep(dd4hep::DetElement& detElement,
                                std::vector<dd4hep::DetElement>& subdetectors,
                                LoggerWrapper logger) {
  const dd4hep::DetElement::Children& children = detElement.children();
  for (auto& child : children) {
    dd4hep::DetElement childDetElement = child.second;
    dd4hep::DetType type{childDetElement.typeFlag()};
    if (childDetElement.type() == "compound") {
      subdetectors.push_back(childDetElement);
      continue;
    }

    if (type.is(dd4hep::DetType::TRACKER)) {
      subdetectors.push_back(childDetElement);
    }
    collectSubDetectors_dd4hep(childDetElement, subdetectors, logger);
  }
}

void collectLayers_dd4hep(dd4hep::DetElement& detElement,
                          std::vector<dd4hep::DetElement>& layers,
                          LoggerWrapper logger) {
  const dd4hep::DetElement::Children& children = detElement.children();
  for (auto& child : children) {
    std::string _expr{"$^"};  // nothing

    dd4hep::rec::VariantParameters* params =
        detElement.extension<dd4hep::rec::VariantParameters>(false);

    if (params != nullptr) {
      _expr = params->value_or<std::string>("layer_pattern", _expr);
      ACTS_VERBOSE("--> Layer pattern for elt " << detElement.name() << ": "
                                                << _expr);
    }
    std::regex expr{_expr};

    dd4hep::DetElement childDetElement = child.second;

    if (std::regex_search(childDetElement.name(), expr)) {
      ACTS_VERBOSE("--> Layer candidate match: " << _expr << " -> "
                                                 << childDetElement.name());
      layers.push_back(childDetElement);
      continue;
    }

    collectLayers_dd4hep(childDetElement, layers, logger);
  }
}

}  // End of namespace Acts

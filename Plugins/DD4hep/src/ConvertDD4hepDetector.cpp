// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/ConvertDD4hepDetector.hpp"
#include <list>
#include <stdexcept>
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/PassiveLayerBuilder.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Plugins/DD4hep/DD4hepLayerBuilder.hpp"
#include "Acts/Plugins/DD4hep/IActsExtension.hpp"
#include "TGeoManager.h"

namespace Acts {
std::unique_ptr<const TrackingGeometry> convertDD4hepDetector(
    dd4hep::DetElement worldDetElement, Logging::Level loggingLevel,
    BinningType bTypePhi, BinningType bTypeR, BinningType bTypeZ,
    double layerEnvelopeR, double layerEnvelopeZ, double defaultLayerThickness,
    const std::function<void(std::vector<dd4hep::DetElement>& detectors)>&
        sortSubDetectors,
    const Acts::GeometryContext& gctx) {
  // create local logger for conversion
  auto DD4hepConverterlogger =
      Acts::getDefaultLogger("DD4hepConversion", loggingLevel);
  ACTS_LOCAL_LOGGER(DD4hepConverterlogger);

  ACTS_INFO("Translating DD4hep geometry into Acts geometry");
  // get the sub detectors of the world detector e.g. beampipe, pixel detector,
  // strip detector
  std::vector<dd4hep::DetElement> subDetectors;
  // go through the detector hierarchies
  collectSubDetectors_dd4hep(worldDetElement, subDetectors);
  // sort to build detector from bottom to top
  sortSubDetectors(subDetectors);
  // the volume builders of the subdetectors
  std::list<std::shared_ptr<const ITrackingVolumeBuilder>> volumeBuilders;
  // the beam pipe volume builder needs special treatment and needs to be added
  // in the end (beampipe exceeds length of all other subdetectors)
  std::shared_ptr<const CylinderVolumeBuilder> beamPipeVolumeBuilder;
  // loop over the sub detectors
  for (auto& subDetector : subDetectors) {
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
  tgbConfig.trackingVolumeBuilders = std::move(volumeFactories);
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
  auto volumeHelper = cylinderVolumeHelper_dd4hep();
  // create local logger for conversion
  auto DD4hepConverterlogger =
      Acts::getDefaultLogger("DD4hepConversion", loggingLevel);
  ACTS_LOCAL_LOGGER(DD4hepConverterlogger);

  Acts::IActsExtension* subDetExtension = nullptr;
  // at this stage not every DetElement needs to have an Extension attached
  try {
    subDetExtension = subDetector.extension<Acts::IActsExtension>();
  } catch (std::runtime_error& e) {
  }
  if (subDetector.type() == "compound") {
    ACTS_VERBOSE("[D] Subdetector : '"
                 << subDetector.name()
                 << "' has no ActsExtension and has type compound - "
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
      ACTS_VERBOSE("[V] Volume : '"
                   << subDetector.name()
                   << "'is a compound volume -> resolve now the sub volumes");

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
      // check if it has a volume extension telling if it is a barrel or an
      // endcap
      IActsExtension* volumeExtension = nullptr;
      try {
        volumeExtension = volumeDetElement.extension<IActsExtension>();
      } catch (std::runtime_error& e) {
        throw std::logic_error(
            std::string("[V] Current DetElement: ") + volumeDetElement.name() +
            std::string(" has no ActsExtension! At this stage it should be a "
                        "detector volume declared as Barrel or Endcap. Please"
                        "check your detector construction."));
      }

      if (volumeExtension->isEndcap()) {
        ACTS_VERBOSE(
            std::string("[V] Subvolume : '") + volumeDetElement.name() +
            std::string("' is a disc volume -> handling as an endcap"));
        if (zPos < 0.) {
          if (nEndCap) {
            throw std::logic_error(
                "[V] Negative Endcap was already given for this "
                "hierachy! Please create a new "
                "DD4hep_SubDetectorAssembly for the next "
                "hierarchy.");
          }
          nEndCap = true;
          ACTS_VERBOSE("[V]       ->is negative endcap");
          collectLayers_dd4hep(volumeDetElement, negativeLayers);
        } else {
          if (pEndCap) {
            throw std::logic_error(
                "[V] Positive Endcap was already given for this "
                "hierachy! Please create a new "
                "DD4hep_SubDetectorAssembly for the next "
                "hierarchy.");
          }
          pEndCap = true;
          ACTS_VERBOSE("[V]       ->is positive endcap");
          collectLayers_dd4hep(volumeDetElement, positiveLayers);
        }
      } else if (volumeExtension->isBarrel()) {
        if (barrel) {
          throw std::logic_error(
              "[V] Barrel was already given for this "
              "hierachy! Please create a new "
              "DD4hep_SubDetectorAssembly for the next "
              "hierarchy.");
        }
        barrel = true;
        ACTS_VERBOSE("[V] Subvolume : "
                     << volumeDetElement.name()
                     << " is a cylinder volume -> handling as a barrel");
        collectLayers_dd4hep(volumeDetElement, centralLayers);
      } else {
        throw std::logic_error(
            std::string("[V] Current DetElement: ") + volumeDetElement.name() +
            std::string(
                " has wrong ActsExtension! At this stage it should be a "
                "detector volume declared as Barrel or Endcap. Please "
                "check your detector construction."));
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
            Acts::getDefaultLogger("SurfaceArrayCreator", loggingLevel));
    // configure LayerCreator
    Acts::LayerCreator::Config lcConfig;
    lcConfig.surfaceArrayCreator = surfaceArrayCreator;
    auto layerCreator = std::make_shared<const Acts::LayerCreator>(
        lcConfig, Acts::getDefaultLogger("LayerCreator", loggingLevel));
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
        lbConfig, Acts::getDefaultLogger("DD4hepLayerBuilder", loggingLevel));

    // get the possible material of the surounding volume
    dd4hep::Material ddmaterial = subDetector.volume().material();
    auto volumeMaterial =
        std::make_shared<const Acts::HomogeneousVolumeMaterial>(Acts::Material(
            ddmaterial.radLength(), ddmaterial.intLength(), ddmaterial.A(),
            ddmaterial.Z(), ddmaterial.density()));

    // the configuration object of the volume builder
    Acts::CylinderVolumeBuilder::Config cvbConfig;

    // Create the sub volume
    // Dimensions are created automatically by adding a tolerance to the
    // layer setup
    cvbConfig.layerEnvelopeR = std::make_pair(layerEnvelopeR, layerEnvelopeR);
    cvbConfig.layerEnvelopeZ = layerEnvelopeZ;
    cvbConfig.trackingVolumeHelper = volumeHelper;
    cvbConfig.volumeSignature = 0;
    cvbConfig.volumeName = subDetector.name();
    cvbConfig.volumeMaterial = volumeMaterial;
    cvbConfig.layerBuilder = dd4hepLayerBuilder;
    auto cylinderVolumeBuilder =
        std::make_shared<const Acts::CylinderVolumeBuilder>(
            cvbConfig,
            Acts::getDefaultLogger("CylinderVolumeBuilder", loggingLevel));
    return cylinderVolumeBuilder;

  } else if ((subDetExtension != nullptr) &&
             (subDetExtension->isPassiveCylinder() ||
              subDetExtension->isBeampipe())) {
    ACTS_VERBOSE("[D] Subdetector : " << subDetector.name()
                                      << " - building a passive cylinder.");
    if (subDetExtension->isBeampipe()) {
      ACTS_VERBOSE("This is the beam pipe - will be built to r->0.");
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
        "[V] Extracting cylindrical volume bounds ( rmin / rmax / "
        "halfZ )=  ( "
        << rMin << " / " << rMax << " / " << halfZ << " )");

    // get the possible material of the surounding volume
    dd4hep::Material ddmaterial = subDetector.volume().material();
    Acts::MaterialProperties pcMaterial(
        ddmaterial.radLength() * UnitConstants::cm,
        ddmaterial.intLength() * UnitConstants::cm, ddmaterial.A(),
        ddmaterial.Z(), ddmaterial.density() / pow(Acts::UnitConstants::cm, 3),
        fabs(tube->GetRmax() - tube->GetRmin()) * UnitConstants::cm);

    // configure the beam pipe layer builder
    Acts::PassiveLayerBuilder::Config bplConfig;
    bplConfig.layerIdentification = subDetector.name();
    bplConfig.centralLayerRadii = std::vector<double>(1, 0.5 * (rMax + rMin));
    bplConfig.centralLayerHalflengthZ = std::vector<double>(1, halfZ);
    bplConfig.centralLayerThickness = std::vector<double>(1, fabs(rMax - rMin));
    bplConfig.centralLayerMaterial = {
        std::make_shared<const HomogeneousSurfaceMaterial>(pcMaterial)};
    auto pcLayerBuilder = std::make_shared<const Acts::PassiveLayerBuilder>(
        bplConfig, Acts::getDefaultLogger(subDetector.name(), loggingLevel));

    // the configuration object of the volume builder
    Acts::CylinderVolumeBuilder::Config cvbConfig;
    cvbConfig.trackingVolumeHelper = volumeHelper;
    cvbConfig.volumeSignature = 0;
    cvbConfig.volumeName = subDetector.name();
    cvbConfig.layerBuilder = pcLayerBuilder;
    cvbConfig.layerEnvelopeR = {layerEnvelopeR, layerEnvelopeR};
    cvbConfig.layerEnvelopeZ = layerEnvelopeZ;
    cvbConfig.buildToRadiusZero = subDetExtension->isBeampipe();

    // beam pipe / passive cylinder volume builder
    auto pcVolumeBuilder = std::make_shared<const Acts::CylinderVolumeBuilder>(
        cvbConfig,
        Acts::getDefaultLogger(
            subDetector.name() + std::string("VolumdeBuilder"), loggingLevel));
    return pcVolumeBuilder;

  } else if ((subDetExtension != nullptr) && subDetExtension->isBarrel()) {
    ACTS_VERBOSE("[D] Subdetector: "
                 << subDetector.name()
                 << " is a (sensitive) Barrel volume - building barrel.");
    /// the dd4hep::DetElements of the layers of the central volume
    std::vector<dd4hep::DetElement> centralLayers;
    collectLayers_dd4hep(subDetector, centralLayers);

    // configure SurfaceArrayCreator
    auto surfaceArrayCreator =
        std::make_shared<const Acts::SurfaceArrayCreator>(
            Acts::getDefaultLogger("SurfaceArrayCreator", loggingLevel));
    // configure LayerCreator
    Acts::LayerCreator::Config lcConfig;
    lcConfig.surfaceArrayCreator = surfaceArrayCreator;
    auto layerCreator = std::make_shared<const Acts::LayerCreator>(
        lcConfig, Acts::getDefaultLogger("LayerCreator", loggingLevel));
    // configure DD4hepLayerBuilder
    Acts::DD4hepLayerBuilder::Config lbConfig;
    lbConfig.configurationName = subDetector.name();
    lbConfig.layerCreator = layerCreator;
    lbConfig.centralLayers = centralLayers;
    lbConfig.bTypePhi = bTypePhi;
    lbConfig.bTypeZ = bTypeZ;
    lbConfig.defaultThickness = defaultLayerThickness;
    auto dd4hepLayerBuilder = std::make_shared<const Acts::DD4hepLayerBuilder>(
        lbConfig, Acts::getDefaultLogger("DD4hepLayerBuilder", loggingLevel));

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

    // get the possible material
    /// @todo volume material currently not used
    dd4hep::Material ddmaterial = subDetector.volume().material();
    auto volumeMaterial =
        std::make_shared<const Acts::HomogeneousVolumeMaterial>(Acts::Material(
            ddmaterial.radLength(), ddmaterial.intLength(), ddmaterial.A(),
            ddmaterial.Z(), ddmaterial.density()));
    cvbConfig.layerEnvelopeR = std::make_pair(layerEnvelopeR, layerEnvelopeR);
    cvbConfig.layerEnvelopeZ = layerEnvelopeZ;
    cvbConfig.trackingVolumeHelper = volumeHelper;
    cvbConfig.volumeSignature = 0;
    cvbConfig.volumeName = subDetector.name();
    cvbConfig.volumeMaterial = volumeMaterial;
    cvbConfig.layerBuilder = dd4hepLayerBuilder;
    auto cylinderVolumeBuilder =
        std::make_shared<const Acts::CylinderVolumeBuilder>(
            cvbConfig,
            Acts::getDefaultLogger("CylinderVolumeBuilder", loggingLevel));
    return cylinderVolumeBuilder;

  } else {
    ACTS_INFO(
        "[D] Subdetector with name : '"
        << subDetector.name()
        << "' has wrong ActsExtension for translation and is not of type "
           "'compound'. If you want to have this DetElement be translated "
           "into the tracking geometry you need add the right "
           "ActsExtension (at this stage the subvolume needs to be "
           "declared as beampipe or barrel) or if it is a compound "
           "DetElement (containing a barrel-endcap hierarchy), the type "
           "needs to be set to 'compound'.");
    return nullptr;
  }
}

std::shared_ptr<const Acts::CylinderVolumeHelper> cylinderVolumeHelper_dd4hep(
    Logging::Level loggingLevel) {
  // create cylindervolumehelper which can be used by all instances
  // hand over LayerArrayCreator
  Acts::LayerArrayCreator::Config lacConfig;
  auto layerArrayCreator = std::make_shared<const Acts::LayerArrayCreator>(
      lacConfig, Acts::getDefaultLogger("LayArrayCreator", loggingLevel));
  // tracking volume array creator
  Acts::TrackingVolumeArrayCreator::Config tvacConfig;
  auto trackingVolumeArrayCreator =
      std::make_shared<const Acts::TrackingVolumeArrayCreator>(
          tvacConfig,
          Acts::getDefaultLogger("TrkVolArrayCreator", loggingLevel));
  // configure the cylinder volume helper
  Acts::CylinderVolumeHelper::Config cvhConfig;
  cvhConfig.layerArrayCreator = layerArrayCreator;
  cvhConfig.trackingVolumeArrayCreator = trackingVolumeArrayCreator;
  auto cylinderVolumeHelper =
      std::make_shared<const Acts::CylinderVolumeHelper>(
          cvhConfig, Acts::getDefaultLogger("CylVolHelper", loggingLevel));

  return cylinderVolumeHelper;
}

void collectCompounds_dd4hep(dd4hep::DetElement& detElement,
                             std::vector<dd4hep::DetElement>& compounds) {
  const dd4hep::DetElement::Children& children = detElement.children();
  for (auto& child : children) {
    dd4hep::DetElement childDetElement = child.second;
    Acts::IActsExtension* detExtension = nullptr;
    try {
      detExtension = childDetElement.extension<Acts::IActsExtension>();
    } catch (std::runtime_error& e) {
    }
    if ((detExtension != nullptr) &&
        (detExtension->isBarrel() || detExtension->isEndcap())) {
      compounds.push_back(childDetElement);
      continue;
    }
    collectCompounds_dd4hep(childDetElement, compounds);
  }
}

void collectSubDetectors_dd4hep(dd4hep::DetElement& detElement,
                                std::vector<dd4hep::DetElement>& subdetectors) {
  const dd4hep::DetElement::Children& children = detElement.children();
  for (auto& child : children) {
    dd4hep::DetElement childDetElement = child.second;
    Acts::IActsExtension* detExtension = nullptr;
    try {
      detExtension = childDetElement.extension<Acts::IActsExtension>();
    } catch (std::runtime_error& e) {
      if (childDetElement.type() == "compound") {
        subdetectors.push_back(childDetElement);
        continue;
      }
    }
    if ((detExtension != nullptr) &&
        (detExtension->isBarrel() || detExtension->isBeampipe())) {
      subdetectors.push_back(childDetElement);
      continue;
    }
    collectSubDetectors_dd4hep(childDetElement, subdetectors);
  }
}

void collectLayers_dd4hep(dd4hep::DetElement& detElement,
                          std::vector<dd4hep::DetElement>& layers) {
  const dd4hep::DetElement::Children& children = detElement.children();
  for (auto& child : children) {
    dd4hep::DetElement childDetElement = child.second;
    Acts::IActsExtension* detExtension = nullptr;
    try {
      detExtension = childDetElement.extension<Acts::IActsExtension>();
    } catch (std::runtime_error& e) {
    }
    if ((detExtension != nullptr) && detExtension->isLayer()) {
      layers.push_back(childDetElement);
      continue;
    }
    collectLayers_dd4hep(childDetElement, layers);
  }
}
}  // End of namespace Acts

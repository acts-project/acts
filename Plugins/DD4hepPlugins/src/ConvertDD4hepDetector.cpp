// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Plugins/DD4hepPlugins/ConvertDD4hepDetector.hpp"
#include <list>
#include <stdexcept>
#include "ACTS/Plugins/DD4hepPlugins/DD4hepLayerBuilder.hpp"
#include "ACTS/Plugins/DD4hepPlugins/IActsExtension.hpp"
#include "ACTS/Tools/CylinderVolumeBuilder.hpp"
#include "ACTS/Tools/CylinderVolumeHelper.hpp"
#include "ACTS/Tools/LayerArrayCreator.hpp"
#include "ACTS/Tools/LayerCreator.hpp"
#include "ACTS/Tools/SurfaceArrayCreator.hpp"
#include "ACTS/Tools/TrackingGeometryBuilder.hpp"
#include "ACTS/Tools/TrackingVolumeArrayCreator.hpp"
#include "TGeoManager.h"

namespace Acts {
std::unique_ptr<TrackingGeometry>
convertDD4hepDetector(DD4hep::Geometry::DetElement worldDetElement,
                      Logging::Level               loggingLevel,
                      BinningType                  bTypePhi,
                      BinningType                  bTypeR,
                      BinningType                  bTypeZ,
                      std::pair<double, double> layerEnvelopeR,
                      double layerEnvelopeZ)
{
  // create local logger for conversion
  auto DD4hepConverterlogger
      = Acts::getDefaultLogger("DD4hepConversion", loggingLevel);
  ACTS_LOCAL_LOGGER(DD4hepConverterlogger);

  ACTS_INFO("Translating DD4hep geometry into ACTS geometry");
  // the return geometry -- and the highest volume
  std::unique_ptr<Acts::TrackingGeometry> trackingGeometry = nullptr;
  // create cylindervolumehelper which can be used by all instances
  // hand over LayerArrayCreator
  auto layerArrayCreator = std::make_shared<Acts::LayerArrayCreator>(
      Acts::getDefaultLogger("LayArrayCreator", loggingLevel));
  // tracking volume array creator
  auto trackingVolumeArrayCreator
      = std::make_shared<Acts::TrackingVolumeArrayCreator>(
          Acts::getDefaultLogger("TrkVolArrayCreator", loggingLevel));
  // configure the cylinder volume helper
  Acts::CylinderVolumeHelper::Config cvhConfig;
  cvhConfig.layerArrayCreator          = layerArrayCreator;
  cvhConfig.trackingVolumeArrayCreator = trackingVolumeArrayCreator;
  auto cylinderVolumeHelper = std::make_shared<Acts::CylinderVolumeHelper>(
      cvhConfig, Acts::getDefaultLogger("CylVolHelper", loggingLevel));

  // get the sub detectors of the world detector e.g. beampipe, pixel detector,
  // strip detector
  std::vector<DD4hep::Geometry::DetElement>     subDetectors;
  const DD4hep::Geometry::DetElement::Children& worldDetectorChildren
      = worldDetElement.children();
  // go through the detector hierarchies
  for (auto& worldDetectorChild : worldDetectorChildren)
    subDetectors.push_back(worldDetectorChild.second);
  // sort by id to build detector from bottom to top
  sort(subDetectors.begin(),
       subDetectors.end(),
       [](const DD4hep::Geometry::DetElement& a,
          const DD4hep::Geometry::DetElement& b) { return (a.id() < b.id()); });

  // the volume builders of the subdetectors
  std::list<std::shared_ptr<ITrackingVolumeBuilder>> volumeBuilders;
  // loop over the sub detectors
  for (auto& subDetector : subDetectors) {
    Acts::IActsExtension* subDetExtension = nullptr;
    // at this stage not every DetElement needs to have an Extension attached
    try {
      subDetExtension = subDetector.extension<Acts::IActsExtension>();
    } catch (std::runtime_error& e) {
      if (subDetector.type() == "compound") {
        ACTS_VERBOSE("[D] Subdetector : '"
                     << subDetector.name()
                     << "' has no ActsExtension and has type compound - "
                        "handling as a compound volume (a hierachy of a "
                        "barrel-endcap structure) and resolving the "
                        "subvolumes...");
        // Now create the Layerbuilders and Volumebuilder
        // the layers
        /// the DD4hep::DetElements of the layers of the negative volume
        std::vector<DD4hep::Geometry::DetElement> negativeLayers;
        /// the DD4hep::DetElements of the layers of the central volume
        std::vector<DD4hep::Geometry::DetElement> centralLayers;
        /// the DD4hep::DetElements of the layers of the positive volume
        std::vector<DD4hep::Geometry::DetElement> positiveLayers;

        // go through sub volumes
        const DD4hep::Geometry::DetElement::Children& subDetectorChildren
            = subDetector.children();
        for (auto& subDetectorChild : subDetectorChildren) {
          DD4hep::Geometry::DetElement volumeDetElement
              = subDetectorChild.second;
          ACTS_VERBOSE(
              "[V] Volume : '"
              << volumeDetElement.name()
              << "'is a compound volume -> resolve now the sub volumes");

          IActsExtension* volumeExtension = nullptr;
          try {
            volumeExtension = volumeDetElement.extension<IActsExtension>();
          } catch (std::runtime_error& e) {
            ACTS_ERROR("[V] Current DetElement: "
                       << volumeDetElement.name()
                       << " has no ActsExtension! At this stage it should be a "
                          "detector volume declared as Barrel or Endcap. Please"
                          "check your detector construction.");
            return nullptr;
          }
          if (volumeExtension->isEndcap()) {
            ACTS_VERBOSE("[V] Subvolume : '"
                         << volumeDetElement.name()
                         << "' is a disc volume -> handling as an endcap");

            if (volumeDetElement.placement()
                    .ptr()
                    ->GetMatrix()
                    ->GetTranslation()[2]
                < 0.) {
              ACTS_VERBOSE("[V]       ->is negative endcap");
              collectLayers(volumeDetElement, negativeLayers);
            } else {
              ACTS_VERBOSE("[V]       ->is positive endcap");
              collectLayers(volumeDetElement, positiveLayers);
            }
          } else if (volumeExtension->isBarrel()) {
            ACTS_VERBOSE("[V] Subvolume : "
                         << volumeDetElement.name()
                         << " is a cylinder volume -> handling as a barrel");
            collectLayers(volumeDetElement, centralLayers);
          } else {
            ACTS_ERROR(
                "[V] Current DetElement: "
                << volumeDetElement.name()
                << " has wrong ActsExtension! At this stage it should be a "
                   "detector volume declared as Barrel or Endcap. Please "
                   "check your detector construction.");
            return nullptr;
          }
        }

        // configure SurfaceArrayCreator
        auto surfaceArrayCreator = std::make_shared<Acts::SurfaceArrayCreator>(
            Acts::getDefaultLogger("SurfaceArrayCreator", loggingLevel));
        // configure LayerCreator
        Acts::LayerCreator::Config lcConfig;
        lcConfig.surfaceArrayCreator = surfaceArrayCreator;
        auto layerCreator            = std::make_shared<Acts::LayerCreator>(
            lcConfig, Acts::getDefaultLogger("LayerCreator", loggingLevel));
        // configure DD4hepLayerBuilder
        Acts::DD4hepLayerBuilder::Config lbConfig;
        lbConfig.configurationName = subDetector.name();
        lbConfig.layerCreator      = layerCreator;
        lbConfig.negativeLayers    = negativeLayers;
        lbConfig.centralLayers     = centralLayers;
        lbConfig.positiveLayers    = positiveLayers;
        lbConfig.bTypePhi          = bTypePhi;
        lbConfig.bTypeR            = bTypeR;
        lbConfig.bTypeZ            = bTypeZ;
        auto dd4hepLayerBuilder    = std::make_shared<Acts::DD4hepLayerBuilder>(
            lbConfig,
            Acts::getDefaultLogger("DD4hepLayerBuilder", loggingLevel));

        // get the possible material of the surounding volume
        DD4hep::Geometry::Material ddmaterial = subDetector.volume().material();
        auto volumeMaterial = std::make_shared<Material>(ddmaterial.radLength(),
                                                         ddmaterial.intLength(),
                                                         ddmaterial.A(),
                                                         ddmaterial.Z(),
                                                         ddmaterial.density());
        // the configuration object of the volume builder
        Acts::CylinderVolumeBuilder::Config cvbConfig;
        cvbConfig.trackingVolumeHelper = cylinderVolumeHelper;
        cvbConfig.volumeSignature      = 0;
        cvbConfig.volumeName           = subDetector.name();
        cvbConfig.volumeMaterial       = volumeMaterial;
        cvbConfig.layerBuilder         = dd4hepLayerBuilder;
        cvbConfig.layerEnvelopeR       = layerEnvelopeR;
        cvbConfig.layerEnvelopeZ       = layerEnvelopeZ;
        auto cylinderVolumeBuilder
            = std::make_shared<Acts::CylinderVolumeBuilder>(
                cvbConfig,
                Acts::getDefaultLogger("CylinderVolumeBuilder", loggingLevel));
        volumeBuilders.push_back(cylinderVolumeBuilder);

      } else {
        ACTS_INFO(
            "[D] Subdetector with name : '"
            << subDetector.name()
            << "' has no ActsExtension and is not of type 'compound'. If you "
               "want to have this DetElement be translated into the tracking "
               "geometry you need add the right ActsExtension (at this stage "
               "the subvolume needs to be declared as beampipe or barrel) or "
               "if it is a compound DetElement (containing a barrel-endcap "
               "hierarchy), the type needs to be set 'compound'.");
        continue;
      }
    }
    if (subDetExtension && subDetExtension->isBeampipe()) {
      ACTS_VERBOSE("[D] Subdetector : "
                   << subDetector.name()
                   << " is the beampipe - building beam pipe.");
      // get the dimensions of the volume
      TGeoShape* geoShape
          = subDetector.placement().ptr()->GetVolume()->GetShape();
      TGeoConeSeg* tube = dynamic_cast<TGeoConeSeg*>(geoShape);
      if (!tube) {
        ACTS_ERROR("Beampipe has wrong shape - needs to be TGeoConeSeg!");
        return nullptr;
      }
      // get the dimension of TGeo and convert lengths
      double rMin  = tube->GetRmin1() * units::_cm;
      double rMax  = tube->GetRmax1() * units::_cm;
      double halfZ = tube->GetDz() * units::_cm;
      ACTS_VERBOSE("[V] Extracting cylindrical volume bounds ( rmin / rmax / "
                   "halfZ )=  ( "
                   << rMin
                   << " / "
                   << rMax
                   << " / "
                   << halfZ
                   << " )");

      // get the possible material of the surounding volume
      DD4hep::Geometry::Material ddmaterial = subDetector.volume().material();
      auto volumeMaterial = std::make_shared<Material>(ddmaterial.radLength(),
                                                       ddmaterial.intLength(),
                                                       ddmaterial.A(),
                                                       ddmaterial.Z(),
                                                       ddmaterial.density());
      // the configuration object of the volume builder
      Acts::CylinderVolumeBuilder::Config cvbConfig;
      cvbConfig.trackingVolumeHelper = cylinderVolumeHelper;
      cvbConfig.volumeSignature      = 0;
      cvbConfig.volumeName           = subDetector.name();
      cvbConfig.buildToRadiusZero    = true;
      cvbConfig.volumeMaterial       = volumeMaterial;
      cvbConfig.volumeDimension      = {rMin, rMax, -halfZ, halfZ};
      auto cylinderVolumeBuilder
          = std::make_shared<Acts::CylinderVolumeBuilder>(
              cvbConfig,
              Acts::getDefaultLogger("CylinderVolumeBuilder", loggingLevel));
      volumeBuilders.push_back(cylinderVolumeBuilder);

    } else if (subDetExtension && subDetExtension->isBarrel()) {
      ACTS_VERBOSE("[D] Subdetector: "
                   << subDetector.name()
                   << " is a Barrel volume - building barrel.");
      /// the DD4hep::DetElements of the layers of the central volume
      std::vector<DD4hep::Geometry::DetElement> centralLayers;
      collectLayers(subDetector, centralLayers);

      // configure SurfaceArrayCreator
      auto surfaceArrayCreator = std::make_shared<Acts::SurfaceArrayCreator>(
          Acts::getDefaultLogger("SurfaceArrayCreator", loggingLevel));
      // configure LayerCreator
      Acts::LayerCreator::Config lcConfig;
      lcConfig.surfaceArrayCreator = surfaceArrayCreator;
      auto layerCreator            = std::make_shared<Acts::LayerCreator>(
          lcConfig, Acts::getDefaultLogger("LayerCreator", loggingLevel));
      // configure DD4hepLayerBuilder
      Acts::DD4hepLayerBuilder::Config lbConfig;
      lbConfig.configurationName = subDetector.name();
      lbConfig.layerCreator      = layerCreator;
      lbConfig.centralLayers     = centralLayers;
      lbConfig.bTypePhi          = bTypePhi;
      lbConfig.bTypeZ            = bTypeZ;
      auto dd4hepLayerBuilder    = std::make_shared<Acts::DD4hepLayerBuilder>(
          lbConfig, Acts::getDefaultLogger("DD4hepLayerBuilder", loggingLevel));

      // get the possible material
      DD4hep::Geometry::Material ddmaterial = subDetector.volume().material();
      auto volumeMaterial = std::make_shared<Material>(ddmaterial.radLength(),
                                                       ddmaterial.intLength(),
                                                       ddmaterial.A(),
                                                       ddmaterial.Z(),
                                                       ddmaterial.density());
      // the configuration object of the volume builder
      Acts::CylinderVolumeBuilder::Config cvbConfig;
      cvbConfig.trackingVolumeHelper = cylinderVolumeHelper;
      cvbConfig.volumeSignature      = 0;
      cvbConfig.volumeName           = subDetector.name();
      cvbConfig.volumeMaterial       = volumeMaterial;
      cvbConfig.layerBuilder         = dd4hepLayerBuilder;
      cvbConfig.layerEnvelopeR       = layerEnvelopeR;
      cvbConfig.layerEnvelopeZ       = layerEnvelopeZ;
      auto cylinderVolumeBuilder
          = std::make_shared<Acts::CylinderVolumeBuilder>(
              cvbConfig,
              Acts::getDefaultLogger("CylinderVolumeBuilder", loggingLevel));

    } else {
      if (subDetector.type() != "compound") {
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
        continue;
      }
      continue;
    }
  }
  // hand over the collected volume builders
  Acts::TrackingGeometryBuilder::Config tgbConfig;
  tgbConfig.trackingVolumeHelper   = cylinderVolumeHelper;
  tgbConfig.trackingVolumeBuilders = volumeBuilders;
  auto trackingGeometryBuilder
      = std::make_shared<Acts::TrackingGeometryBuilder>(tgbConfig);
  return (trackingGeometryBuilder->trackingGeometry());
}

void
collectLayers(DD4hep::Geometry::DetElement&              detElement,
              std::vector<DD4hep::Geometry::DetElement>& layers)
{
  const DD4hep::Geometry::DetElement::Children& children
      = detElement.children();
  if (!children.empty()) {
    for (auto& child : children) {
      DD4hep::Geometry::DetElement childDetElement = child.second;
      Acts::IActsExtension*        detExtension    = nullptr;
      try {
        detExtension = childDetElement.extension<Acts::IActsExtension>();
      } catch (std::runtime_error& e) {
        continue;
      }
      if (detExtension->isLayer()) layers.push_back(childDetElement);
      collectLayers(childDetElement, layers);
    }
  }
}
}  // End of namespace Acts

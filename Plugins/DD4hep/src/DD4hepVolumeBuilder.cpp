// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/DD4hep/DD4hepVolumeBuilder.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Root/TGeoPrimitivesHelper.hpp"

#include <stdexcept>
#include <utility>

#include <DD4hep/Alignments.h>
#include <DD4hep/DetElement.h>
#include <DD4hep/Volumes.h>
#include <RtypesCore.h>

using namespace Acts;

namespace ActsPlugins {

DD4hepVolumeBuilder::DD4hepVolumeBuilder(
    const DD4hepVolumeBuilder::Config& config,
    std::unique_ptr<const Logger> logger)
    : m_cfg(), m_logger(std::move(logger)) {
  setConfiguration(config);
}

DD4hepVolumeBuilder::~DD4hepVolumeBuilder() = default;

void DD4hepVolumeBuilder::setConfiguration(
    const DD4hepVolumeBuilder::Config& config) {
  m_cfg = config;
}

std::vector<std::shared_ptr<TrackingVolume>>
DD4hepVolumeBuilder::centralVolumes() const {
  if (m_cfg.centralVolumes.empty()) {
    ACTS_VERBOSE("[L] No layers handed over for central volume!");
    return {};
  }

  ACTS_VERBOSE(
      "[L] Received layers for central volume -> creating "
      "cylindrical layers");

  // Resulting volumes
  MutableTrackingVolumeVector volumes;
  // Inner/outer radius and half length of the barrel
  double rMin = 0, rMax = 0, dz = 0;

  // Go through volumes
  for (auto& detElement : m_cfg.centralVolumes) {
    // Access the global transformation matrix of the volume
    auto transform =
        convertTransform(&(detElement.nominal().worldTransformation()));
    // Get the shape of the volume
    TGeoShape* geoShape = detElement.placement().ptr()->GetVolume()->GetShape();

    if (geoShape != nullptr) {
      TGeoTubeSeg* tube = dynamic_cast<TGeoTubeSeg*>(geoShape);
      if (tube == nullptr) {
        ACTS_ERROR(
            "[L] Cylinder layer has wrong shape - needs to be TGeoTubeSeg!");
        throw std::logic_error{
            "[L] Cylinder layer has wrong shape - needs to be TGeoTubeSeg!"};
      }

      // Extract the boundaries
      rMin = tube->GetRmin() * UnitConstants::cm;
      rMax = tube->GetRmax() * UnitConstants::cm;
      dz = tube->GetDz() * UnitConstants::cm;

    } else {
      throw std::logic_error(
          std::string("Volume DetElement: ") + detElement.name() +
          std::string(" has not a shape "
                      "added to its extension. Please check your detector "
                      "constructor!"));
    }
    // Build boundaries
    volumes.push_back(std::make_shared<TrackingVolume>(
        transform, std::make_shared<CylinderVolumeBounds>(rMin, rMax, dz)));
  }
  return volumes;
}

Transform3 DD4hepVolumeBuilder::convertTransform(
    const TGeoMatrix* tGeoTrans) const {
  // Get the placement and orientation in respect to its mother
  const Double_t* rotation = tGeoTrans->GetRotationMatrix();
  const Double_t* translation = tGeoTrans->GetTranslation();
  return TGeoPrimitivesHelper::makeTransform(
      Vector3(rotation[0], rotation[3], rotation[6]),
      Vector3(rotation[1], rotation[4], rotation[7]),
      Vector3(rotation[2], rotation[5], rotation[8]),
      Vector3(translation[0] * UnitConstants::cm,
              translation[1] * UnitConstants::cm,
              translation[2] * UnitConstants::cm));
}

}  // namespace ActsPlugins

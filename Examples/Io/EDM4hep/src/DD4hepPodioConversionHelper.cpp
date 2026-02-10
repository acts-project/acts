// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/DD4hepPodioConversionHelper.hpp"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsPlugins/DD4hep/DD4hepDetectorElement.hpp"
#include "ActsPodioEdm/TrackerHitLocalCollection.h"

#include <DD4hep/Detector.h>
#include <DD4hep/Volumes.h>

namespace ActsExamples {

DD4hepPodioConversionHelper::DD4hepPodioConversionHelper(
    const DD4hepDetector& detector,
    const ActsPodioEdm::TrackerHitLocalCollection& trackerHitLocalCollection)
    : m_detector(&detector),
      m_trackerHitLocalCollection(&trackerHitLocalCollection) {}

std::optional<ActsPlugins::PodioUtil::Identifier>
DD4hepPodioConversionHelper::surfaceToIdentifier(
    const Acts::Surface& surface) const {
  const auto* detElement =
      dynamic_cast<const ActsPlugins::DD4hepDetectorElement*>(
          surface.surfacePlacement());
  if (detElement == nullptr) {
    return std::nullopt;
  }

  return static_cast<ActsPlugins::PodioUtil::Identifier>(
      detElement->sourceElement().volumeID());
}

const Acts::Surface* DD4hepPodioConversionHelper::identifierToSurface(
    ActsPlugins::PodioUtil::Identifier identifier) const {
  auto detElement =
      m_detector->dd4hepDetector().volumeManager().lookupDetElement(identifier);
  if (!detElement.isValid()) {
    return nullptr;
  }

  auto* dd4hepDetElementExtension =
      detElement.extension<ActsPlugins::DD4hepDetectorElementExtension>();
  if (dd4hepDetElementExtension == nullptr) {
    return nullptr;
  }

  return &dd4hepDetElementExtension->detectorElement().surface();
}

ActsPodioEdm::TrackerHitLocal
DD4hepPodioConversionHelper::sourceLinkToTrackerHitLocal(
    const Acts::SourceLink& sourceLink) const {
  const auto* indexSourceLink = sourceLink.getPtr<IndexSourceLink>();
  if (indexSourceLink == nullptr) {
    throw std::invalid_argument("Source link is not an IndexSourceLink");
  }
  // This assumes that the conversion was index preserving!
  return m_trackerHitLocalCollection->at(indexSourceLink->index());
}

}  // namespace ActsExamples

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DetectorCommons/StructureSelector.hpp"

#include "ActsExamples/Utilities/GroupBy.hpp"

namespace {

struct SensitiveGetter {
  std::vector<std::shared_ptr<const Acts::Surface>> selected;
  /// @brief Visitor call operator
  /// @param surface
  void operator()(const Acts::Surface* surface) {
    if (surface != nullptr) {
      auto geoId = surface->geometryId();
      if (geoId.sensitive() != 0u ||
          surface->associatedDetectorElement() != nullptr) {
        selected.emplace_back(surface->getSharedPtr());
      }
    }
  }
};

}  // namespace

ActsExamples::StructureSelector::StructureSelector(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry)
    : m_trackingGeometry(std::move(trackingGeometry)) {
  if (!m_trackingGeometry) {
    throw std::runtime_error("Tracking geometry is not provided");
  }
  SensitiveGetter getter;
  m_trackingGeometry->visitSurfaces(getter);
  m_surfaceMultiSet = GeometryIdMultiset<std::shared_ptr<const Acts::Surface>>(
      getter.selected.begin(), getter.selected.end());
}

std::vector<std::shared_ptr<const Acts::Surface>>
ActsExamples::StructureSelector::selectSurfaces(
    const Acts::GeometryIdentifier& geoId) const {
  auto selectedRange =
      selectLowestNonZeroGeometryObject(m_surfaceMultiSet, geoId);
  return {selectedRange.begin(), selectedRange.end()};
}

std::unordered_map<Acts::GeometryIdentifier, Acts::Transform3>
ActsExamples::StructureSelector::selectedTransforms(
    const Acts::GeometryContext& gctx,
    const Acts::GeometryIdentifier& geoId) const {
  std::unordered_map<Acts::GeometryIdentifier, Acts::Transform3> transforms;
  auto selectedSurfaces = selectSurfaces(geoId);
  for (const auto& surface : selectedSurfaces) {
    transforms[surface->geometryId()] = surface->transform(gctx);
  }
  return transforms;
}

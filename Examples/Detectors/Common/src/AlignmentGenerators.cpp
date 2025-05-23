// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DetectorCommons/AlignmentGenerators.hpp"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TransformStore.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"

namespace {

std::unordered_map<Acts::GeometryIdentifier, Acts::Transform3> selectTransforms(
    const Acts::GeometryContext& gctx, const Acts::TrackingGeometry& tGeometry,
    const std::vector<Acts::GeometryIdentifier>& selection) {
  // The selection is not empty, so we need to select the
  // transforms for the selected elements
  // This is a bit of a hack, but it works for now

  // Run through the entire hierarchy map and fill it into a
  // convenient multimap
  std::unordered_map<Acts::GeometryIdentifier, Acts::Transform3>
      identifiedTransforms;

  ActsExamples::GeometryIdMultiset<const Acts::Surface*> surfaceMultiSet;
  for (const auto& [geoId, surface] : tGeometry.geoIdSurfaceMap()) {
    if (selection.empty()) {
      identifiedTransforms.emplace_hint(identifiedTransforms.end(), geoId,
                                        surface->transform(gctx));
    } else {
      surfaceMultiSet.emplace_hint(surfaceMultiSet.end(), surface);
    }
  }

  // The path of the selection
  if (!selection.empty()) {
    for (const auto& sId : selection) {
      auto selectedSurfaces =
          ActsExamples::selectLowestNonZeroGeometryObject(surfaceMultiSet, sId);
      for (const auto& surface : selectedSurfaces) {
        identifiedTransforms.emplace(surface->geometryId(),
                                     surface->transform(gctx));
      }
    }
  }
  return identifiedTransforms;
}

}  // namespace

ActsExamples::GlobalShift::GlobalShift(
    const Acts::GeometryContext& gctx,
    const Acts::TrackingGeometry& trackingGeometry,
    const std::vector<Acts::GeometryIdentifier>& selection,
    const Acts::Vector3& ishift)
    : shift(ishift),
      nominalTransforms(selectTransforms(gctx, trackingGeometry, selection)) {
  // Loop over the selected elements and store the nominal transforms
}

std::shared_ptr<Acts::ITransformStore> ActsExamples::GlobalShift::operator()() {
  // Create the transform store
  auto contextualTransforms = nominalTransforms;
  std::for_each(contextualTransforms.begin(), contextualTransforms.end(),
                [this](auto& itrf) {
                  // Apply the shift to the transform
                  itrf.second = itrf.second * shift;
                });

  return std::make_shared<Acts::TransformStoreGeometryId>(contextualTransforms);
}

ActsExamples::PerpendicularScale::PerpendicularScale(
    const Acts::GeometryContext& gctx,
    const Acts::TrackingGeometry& trackingGeometry,
    const std::vector<Acts::GeometryIdentifier>& selection, double iexpansion)
    : expansion(iexpansion),
      nominalTransforms(selectTransforms(gctx, trackingGeometry, selection)) {
  // Loop over the selected elements and store the nominal transforms
}

std::shared_ptr<Acts::ITransformStore> ActsExamples::PerpendicularScale::
    operator()() {
  // Create the transform store
  auto contextualTransforms = nominalTransforms;
  std::for_each(contextualTransforms.begin(), contextualTransforms.end(),
                [this](auto& itrf) {
                  // Apply the radial expansion to the transform
                  itrf.second.translation()[0] *= expansion;
                  itrf.second.translation()[1] *= expansion;});

  return std::make_shared<Acts::TransformStoreGeometryId>(contextualTransforms);
}

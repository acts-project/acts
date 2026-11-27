// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsExamples/DetectorCommons/AlignmentContext.hpp"

#include <stdexcept>

namespace ActsExamples {

/// Updater that writes each aligned transform into @p store.
///
/// AlignmentAlgorithm calls this from finalize(), after event processing has
/// finished. The store is not synchronized for concurrent updates.
inline ActsAlignment::AlignedTransformUpdater makeAlignedTransformUpdater(
    const std::shared_ptr<GeoIdAlignmentStore>& store) {
  if (store == nullptr) {
    throw std::invalid_argument(
        "makeAlignedTransformUpdater: alignment store is null");
  }

  return ActsAlignment::AlignedTransformUpdater(
      [store](Acts::SurfacePlacementBase* placement,
              const Acts::GeometryContext& /*gctx*/,
              const Acts::Transform3& aTransform) -> bool {
        if (placement == nullptr) {
          return false;
        }

        const Acts::Surface& surface = placement->surface();
        Acts::GeometryIdentifier geoId = surface.geometryId();
        store->setTransform(geoId, aTransform);

        return true;
      });
}
}  // namespace ActsExamples

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsExamples/DetectorCommons/Aligned.hpp"
#include "ActsExamples/DetectorCommons/AlignmentContext.hpp"

namespace ActsExamples {

/// @brief Create an aligned transform updater with an alignment store
/// @param store An alignment store to update with aligned transforms
/// @return An AlignedTransformUpdater that updates the provided store
///
/// Usage example:
/// @code
/// // 1. Create a store (this will be updated during alignment)
/// auto alignmentStore = std::make_shared<GeoIdAlignmentStore>(
///     std::unordered_map<Acts::GeometryIdentifier, Acts::Transform3>{});
///
/// // 2. Create the updater that will modify the store
/// auto updater = makeAlignedTransformUpdater(alignmentStore);
/// alignmentConfig.alignedTransformUpdater = updater;
///
/// // 3. Set up the AlignmentDecorator to use this store
/// // The decorator will provide an AlignmentContext pointing to this store
/// AlignmentDecorator::Config decoratorConfig;
/// // Option A: Use the store directly for all events
/// decoratorConfig.iovStores = {{{0, 999999}, alignmentStore}};
///
/// // 4. During alignment iterations:
/// // - AlignmentAlgorithm reads detector positions via AlignmentContext
/// // - AlignmentAlgorithm updates alignmentStore via alignedTransformUpdater
/// // - Next iteration sees the updated transforms automatically
/// @endcode
inline ActsAlignment::AlignedTransformUpdater makeAlignedTransformUpdater(
    const std::shared_ptr<GeoIdAlignmentStore>& store) {
  return ActsAlignment::AlignedTransformUpdater(
      [store](Acts::SurfacePlacementBase* placement,
              const Acts::GeometryContext& /*gctx*/,
              const Acts::Transform3& aTransform) -> bool {
        if (store == nullptr || placement == nullptr) {
          return false;
        }

        const Acts::Surface& surface = placement->surface();
        Acts::GeometryIdentifier geoId = surface.geometryId();
        store->setTransform(geoId, aTransform);

        return true;
      });
}

/// @brief Create an aligned transform updater that captures transforms to a map
/// @param transformMap A map to store the aligned transforms
/// @return An AlignedTransformUpdater that populates the provided map
///
/// This is useful when you want to collect the alignment results
/// without creating a full alignment store
inline ActsAlignment::AlignedTransformUpdater makeAlignedTransformUpdater(
    const std::shared_ptr<std::unordered_map<Acts::GeometryIdentifier,
                                             Acts::Transform3>>& transformMap) {
  return ActsAlignment::AlignedTransformUpdater(
      [transformMap](Acts::SurfacePlacementBase* placement,
                     const Acts::GeometryContext& /*gctx*/,
                     const Acts::Transform3& aTransform) -> bool {
        if (transformMap == nullptr || placement == nullptr) {
          return false;
        }

        const Acts::Surface& surface = placement->surface();
        Acts::GeometryIdentifier geoId = surface.geometryId();
        (*transformMap)[geoId] = aTransform;

        return true;
      });
}
}  // namespace ActsExamples

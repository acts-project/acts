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

#include <stdexcept>

namespace ActsExamples {

/// @brief Create an aligned transform updater with an alignment store
/// @param store An alignment store to update with aligned transforms
/// @return An AlignedTransformUpdater that updates the provided store
///
/// @warning Single-threaded only. The updater writes the shared store from
/// execute(), and other events read it through contextualTransform. The source
/// cannot force Sequencer to one thread, so the alignment job must be run with
/// numThreads = 1. A later change should make this write safe for more threads.
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

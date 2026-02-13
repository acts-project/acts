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

/// @brief Mutable alignment store that supports updates
/// This extends GeoIdAlignmentStore with the ability to update transforms
class MutableGeoIdAlignmentStore : public GeoIdAlignmentStore {
 public:
  using GeoIdAlignmentStore::GeoIdAlignmentStore;

  /// Update or insert a transform for a given geometry ID
  /// @param geoId The geometry identifier
  /// @param transform The new transform to store
  void setTransform(Acts::GeometryIdentifier geoId,
                    const Acts::Transform3& transform) {
    m_transformMap[geoId] = transform;
  }

  /// Get mutable access to the transform map
  std::unordered_map<Acts::GeometryIdentifier, Acts::Transform3>&
  getTransformMap() {
    return m_transformMap;
  }

  /// Get const access to the transform map
  const std::unordered_map<Acts::GeometryIdentifier, Acts::Transform3>&
  getTransformMap() const {
    return m_transformMap;
  }

 protected:
  // Allow access to base class's private members
  using GeoIdAlignmentStore::m_transformMap;
};

/// @brief Default aligned transform updater
/// This is a no-op updater that just returns true
/// Use makeAlignedTransformUpdater() to create a functional updater
inline ActsAlignment::AlignedTransformUpdater alignedTransformUpdater =
    ActsAlignment::AlignedTransformUpdater(
        [](Acts::SurfacePlacementBase*
           /*placement*/,
           const Acts::GeometryContext& /*gctx*/, const Acts::Transform3&
           /*aTransform*/) -> bool {
          // Default no-op updater
          // Use makeAlignedTransformUpdater() for actual functionality
          return true;
        });

/// @brief Create an aligned transform updater with a mutable store
/// @param store A mutable alignment store to update with aligned transforms
/// @return An AlignedTransformUpdater that updates the provided store
///
/// Usage example:
/// @code
/// // 1. Create a mutable store (this will be updated during alignment)
/// auto alignmentStore = std::make_shared<MutableGeoIdAlignmentStore>(
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
///
/// @note The key insight: Since AlignmentContext holds a pointer to the store,
/// any updates to the store are immediately visible in subsequent iterations.
/// The store object must remain alive throughout the alignment process.
inline ActsAlignment::AlignedTransformUpdater makeAlignedTransformUpdater(
    std::shared_ptr<MutableGeoIdAlignmentStore> store) {
  return ActsAlignment::AlignedTransformUpdater(
      [store](Acts::SurfacePlacementBase* placement,
              const Acts::GeometryContext& /*gctx*/,
              const Acts::Transform3& aTransform) -> bool {
        if (store == nullptr || placement == nullptr) {
          return false;
        }

        // Get the geometry ID from the surface placement's surface
        const Acts::Surface& surface = placement->surface();
        Acts::GeometryIdentifier geoId = surface.geometryId();

        // Update the store with the new aligned transform
        // This update is immediately visible to any code using the store
        // via AlignmentContext, enabling iterative alignment
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
    std::shared_ptr<
        std::unordered_map<Acts::GeometryIdentifier, Acts::Transform3>>
        transformMap) {
  return ActsAlignment::AlignedTransformUpdater(
      [transformMap](Acts::SurfacePlacementBase* placement,
                     const Acts::GeometryContext& /*gctx*/,
                     const Acts::Transform3& aTransform) -> bool {
        if (transformMap == nullptr || placement == nullptr) {
          return false;
        }

        // Get the geometry ID from the surface placement's surface
        const Acts::Surface& surface = placement->surface();
        Acts::GeometryIdentifier geoId = surface.geometryId();

        // Store the aligned transform in the map
        (*transformMap)[geoId] = aTransform;

        return true;
      });
}
}  // namespace ActsExamples



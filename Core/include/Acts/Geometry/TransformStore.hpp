// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <unordered_map>

namespace Acts {

/// @brief Base class for alignment stores which can retrieve contextual
/// transforms from specific surfaces.
///
/// Possible implementations may take the geometry identifier, access the
/// detector element or use other means to identify which transform is being
/// queired.
class ITransformStore {
 public:
  /// @brief  Virtual destructor
  virtual ~ITransformStore() = default;

  /// Retrieve the contextual transform for a given surface
  ///
  /// @param surface the  surface for which the contextual tranfrom is requested
  /// @return a pointer to the transform if found, otherwise nullptr
  virtual const Transform3* contextualTransform(
      const Surface& surface) const = 0;
};

/// One possible implementation with a simple unordered map that relates a
/// surface via its geometry identifier to a transform within the store
///
/// To use this store, the GeometryContext of the corresponding geometry must
/// be decorated with such a store and equipped to use it.
class TransformStoreGeometryId : public ITransformStore {
 public:
  /// Constructor from an unordered map of geometry ids and transforms
  /// @param transformMap the map of geometry ids and transforms
  explicit TransformStoreGeometryId(
      std::unordered_map<GeometryIdentifier, Acts::Transform3> transformMap)
      : m_identifiedTransforms(std::move(transformMap)) {}

  /// @copydoc ITransformStore::contextualTransform
  const Transform3* contextualTransform(const Surface& surface) const override {
    auto it = m_identifiedTransforms.find(surface.geometryId());
    if (it != m_identifiedTransforms.end()) {
      return &(it->second);
    }
    return nullptr;
  }

 private:
  /// The geometry id map
  std::unordered_map<GeometryIdentifier, Transform3> m_identifiedTransforms;
};

}  // namespace Acts

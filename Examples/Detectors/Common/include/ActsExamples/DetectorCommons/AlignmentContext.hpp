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

namespace ActsExamples {

/// @brief Base class for alignment stores which can retrieve contextual
/// transforms from specific surfaces.
///
/// Possible implementations may take the geometry identifier, access the
/// detector element or use other means to identify which transform is being
/// queried.
class IAlignmentStore {
 public:
  /// @brief  Virtual destructor
  virtual ~IAlignmentStore() = default;

  /// Retrieve the contextual transform for a given surface
  ///
  /// @param surface the  surface for which the contextual transform is requested
  /// @return a pointer to the transform if found, otherwise nullptr
  virtual const Acts::Transform3* contextualTransform(
      const Acts::Surface& surface) const = 0;
};

/// A simple struct holding a store raw pointer, ownership should not be in the
/// Context as the store may expand lifetime beyond a context scope
struct AlignmentContext {
  /// The store pointer
  const IAlignmentStore* store{nullptr};
};

/// One possible implementation with a simple unordered map that relates a
/// surface via its geometry identifier to a transform within the store
///
/// To use this store, the GeometryContext of the corresponding geometry must
/// be decorated with such a store and equipped to use it.
class GeoIdAlignmentStore : public IAlignmentStore {
 public:
  /// Constructor from an unordered map of geometry ids and transforms
  /// @param transformMap the map of geometry ids and transforms
  explicit GeoIdAlignmentStore(
      std::unordered_map<Acts::GeometryIdentifier, Acts::Transform3>
          transformMap)
      : m_transformMap(std::move(transformMap)) {}

  /// @copydoc IAlignmentStore::contextualTransform
  const Acts::Transform3* contextualTransform(
      const Acts::Surface& surface) const override {
    auto it = m_transformMap.find(surface.geometryId());
    if (it != m_transformMap.end()) {
      return &(it->second);
    }
    return nullptr;
  }

 private:
  /// The geometry id map
  std::unordered_map<Acts::GeometryIdentifier, Acts::Transform3> m_transformMap;
};

}  // namespace ActsExamples

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <unordered_map>

namespace Acts {

/// @brief  Base class for alignment stores which is valid within one
/// geometry context.
class IDD4hepAlignmentStore {
 public:
  /// @brief  Virtual destructor
  virtual ~IDD4hepAlignmentStore() = default;

  /// @brief  Get the alignment transform for a given DD4hep detector element
  /// @param dd4hepElement the DD4hep detector element for which the alignment is queried
  /// @return a pointer to the transform if found, otherwise nullptr
  virtual const Transform3* contextualTransform(
      const DD4hepDetectorElement& dd4hepElement) const = 0;
};

/// @brief A simple alignment store for the DD4hep geometry with a lookup map
/// from geometry id, that is taken from the DD4hep detector element via the
/// Acts::Surface
class DD4hepAlignmentStoreGeometryId : public IDD4hepAlignmentStore {
 public:
  /// Constructor from a map of geometry ids and transforms
  /// @param transformMap the map of geometry ids and transforms
  DD4hepAlignmentStoreGeometryId(
      std::unordered_map<GeometryIdentifier, Acts::Transform3> transformMap)
      : m_identifiedTransforms(std::move(transformMap)) {}

  /// @brief Get the alignment transform for a given DD4hep detector element
  /// @param dd4hepElement the DD4hep detector element for which the alignment is queried
  /// @return a pointer to the transform if found, otherwise nullptr
  const Transform3* contextualTransform(
      const DD4hepDetectorElement& dd4hepElement) const {
    auto it = m_identifiedTransforms.find(dd4hepElement.surface().geometryId());
    if (it != m_identifiedTransforms.end()) {
      return &(it->second);
    }
    return nullptr;
  }

 private:
  /// The geometry id map
  std::unordered_map<GeometryIdentifier, Acts::Transform3>
      m_identifiedTransforms;
};

}  // namespace Acts

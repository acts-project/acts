// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"

#include <memory>
#include <unordered_map>
#include <vector>

namespace Acts {
class Surface;
}

namespace ActsExamples {

/// StructureSelector is a utility class to select a specific structure
class StructureSelector {
 public:
  /// Constructor
  /// @param trackingGeometry The tracking geometry to select from
  explicit StructureSelector(
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry);

  /// Select a structure by its identifier
  /// @param geoId The identifier of the structure to select
  /// @return A vector of pointers to the selected surfaces
  std::vector<std::shared_ptr<const Acts::Surface>> selectSurfaces(
      const Acts::GeometryIdentifier& geoId) const;

  /// Select a structure by its identified transforms
  /// @param gctx The geometry context
  /// @param geoId The geometry identifier of the structure to select
  /// @return an unordered map of geometry identifiers to transforms
  std::unordered_map<Acts::GeometryIdentifier, Acts::Transform3>
  selectedTransforms(const Acts::GeometryContext& gctx,
                     const Acts::GeometryIdentifier& geoId) const;

 private:
  std::shared_ptr<const Acts::TrackingGeometry> m_trackingGeometry;
  GeometryIdMultiset<std::shared_ptr<const Acts::Surface>> m_surfaceMultiSet;
};
}  // namespace ActsExamples

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
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "ActsTests/CommonHelpers/DetectorElementStub.hpp"
#include "ActsTests/CommonHelpers/PredefinedMaterials.hpp"

#include <functional>
#include <memory>
#include <vector>

namespace ActsTests {

struct CubicTrackingGeometry {
  /// Default constructor for the Cubic tracking geometry
  ///
  /// @param gctx the geometry context for this geometry at building time
  explicit CubicTrackingGeometry(const Acts::GeometryContext& gctx);

  /// Call operator to build the standard cubic tracking geometry
  std::shared_ptr<const Acts::TrackingGeometry> operator()();

  Acts::RotationMatrix3 rotation = Acts::RotationMatrix3::Identity();
  std::shared_ptr<const Acts::RectangleBounds> rBounds = nullptr;
  std::shared_ptr<const Acts::ISurfaceMaterial> surfaceMaterial = nullptr;

  std::vector<std::unique_ptr<const DetectorElementStub>> detectorStore = {};

  std::reference_wrapper<const Acts::GeometryContext> geoContext;
};
}  // namespace ActsTests

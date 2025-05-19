// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CylinderLayer.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/PassiveLayerBuilder.hpp"
#include "Acts/Geometry/ProtoLayer.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"

#include <functional>
#include <numbers>
#include <vector>

namespace Acts::Test {

struct CylindricalTrackingGeometry {
  std::reference_wrapper<const GeometryContext> geoContext;

  /// Only allowed constructor with reference wrapper
  explicit CylindricalTrackingGeometry(const GeometryContext& gctx)
      : geoContext(gctx) {}

  using DetectorStore = std::vector<std::unique_ptr<const DetectorElementStub>>;

  /// The detector store for memory management
  DetectorStore detectorStore = {};

  /// Generator of surfaces for a ring
  ///
  /// @param detStore The DetectorStore for storing the modules
  /// @param moduleHalfXminY The half length in X (at Y min) of the module
  /// @param moduleHalfXmaxY The half length in X (at Y max) of the module
  /// @param moduleHalfY The half length in Y of the module
  /// @param moduleThickness The module thickness
  /// @param moduleTilt The tilt out of the plane for discs
  /// @param ringRadius The central radius of the ring
  /// @param ringZ The z position of the ring
  /// @param zStagger The z offset of phi modules
  /// @param nPhi The number of phi modules
  ///
  /// @return A vector of Surfaces
  std::vector<const Surface*> surfacesRing(
      DetectorStore& detStore, double moduleHalfXminY, double moduleHalfXmaxY,
      double moduleHalfY, double moduleThickness, double moduleTilt,
      double ringRadius, double ringZ, double zStagger, int nPhi);

  /// Generator of surfaces for a cylindrical layer
  ///
  /// @param detStore The DetectorStore for storing the modules
  /// @param moduleHalfX The half length in X of the module
  /// @param moduleHalfY The half length in Y of the module
  /// @param moduleThickness The module thickness
  /// @param moduleTilePhi The tilt in phi direction of the module
  /// @param layerRadius The radius of the cylindrical layer
  /// @param radialStagger The radial delta of modules next in z
  /// @param longitudinalOverlap The z overlap of modules next in z
  /// @param binningSchema The number of bins in phi/z
  ///
  /// @return A vector of Surfaces
  std::vector<const Surface*> surfacesCylinder(
      DetectorStore& detStore, double moduleHalfX, double moduleHalfY,
      double moduleThickness, double moduleTiltPhi, double layerRadius,
      double radialStagger, double longitudinalOverlap,
      const std::pair<int, int>& binningSchema);

  /// Helper method for cylinder layer
  /// create the positions for module surfaces on a cylinder
  std::vector<Vector3> modulePositionsCylinder(
      double radius, double zStagger, double moduleHalfLength, double lOverlap,
      const std::pair<int, int>& binningSchema);

  // @brief Call operator for the creation method of the tracking geometry
  std::shared_ptr<const TrackingGeometry> operator()();
};

}  // namespace Acts::Test

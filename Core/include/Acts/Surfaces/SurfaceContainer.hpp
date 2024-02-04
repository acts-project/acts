// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Detector/Detector.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"

#include <vector>
namespace Acts {

using SurfacePtrsContainer = std::vector<const Acts::Surface*>;
using DetectorPtr = std::shared_ptr<Acts::Experimental::Detector>;
using TrackingGeometryPtr = std::shared_ptr<const Acts::TrackingGeometry>;

/// @brief Class that extracts a vector of Surface pointers from
/// either a detector object OR a tracking geometry object

class SurfaceContainer {
 private:
  SurfacePtrsContainer m_surfacePtrs;

  // create a struct that utilizes surfaceVisitor method of Tracking Geometry
  struct SurfaceVisitor {
    SurfacePtrsContainer surfacePtrs;
    void operator()(const Acts::Surface* surface) {
      surfacePtrs.push_back(surface);
    }
  };

  SurfacePtrsContainer getPtrs(const DetectorPtr& detector) const;

  SurfacePtrsContainer getPtrs(const TrackingGeometryPtr& tGeometryPtr) const;

 public:
  /// Constructor from detector object
  ///
  /// @param detector Shared pointer to detector containing the surfaces
  SurfaceContainer(const DetectorPtr& detector)
      : m_surfacePtrs(getPtrs(detector)) {}

  /// Constructor from tracking geometry object
  ///
  /// @param tGeometryPtr Shared pointer to tracking geometry
  /// containing the surfaces
  SurfaceContainer(const TrackingGeometryPtr& tGeometryPtr)
      : m_surfacePtrs(getPtrs(tGeometryPtr)) {}

  SurfaceContainer(const SurfacePtrsContainer& surfaceVec)
      : m_surfacePtrs(surfaceVec) {}

  // Get surface pointer vector
  SurfacePtrsContainer surfacePtrs() const { return m_surfacePtrs; }
};

}  // namespace Acts

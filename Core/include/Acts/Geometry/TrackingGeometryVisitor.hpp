// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <functional>

namespace Acts {

class TrackingVolume;
class Portal;
class Surface;
class Layer;

template <typename T>
class BoundarySurfaceT;

/// @brief Visitor interface for traversing the tracking geometry hierarchy
///
/// This visitor allows for const access to traverse and inspect the tracking
/// geometry components without modifying them. It's used for operations like
/// visualization, validation, or collecting information about the geometry
/// structure.
class TrackingGeometryVisitor {
 public:
  virtual ~TrackingGeometryVisitor() = 0;

  /// @brief Visit a tracking volume in the geometry
  /// @param volume The tracking volume being visited
  /// @note Called for each volume in the geometry hierarchy during traversal
  virtual void visitVolume(const TrackingVolume& volume);

  /// @brief Visit a portal (boundary between volumes)
  /// @param portal The portal being visited
  /// @note Called for each portal encountered during geometry traversal
  virtual void visitPortal(const Portal& portal);

  /// @brief Visit a surface in the geometry
  /// @param surface The surface being visited
  /// @note Called for each surface encountered during geometry traversal
  virtual void visitSurface(const Surface& surface);

  // Gen 1
  /// @brief Visit a detector layer
  /// @param layer The layer being visited
  /// @note Called for each layer encountered during geometry traversal
  virtual void visitLayer(const Layer& layer);

  /// @brief Visit a boundary surface between tracking volumes
  /// @param boundary The boundary surface being visited
  /// @note Called for each boundary surface encountered during geometry traversal
  virtual void visitBoundarySurface(
      const BoundarySurfaceT<TrackingVolume>& boundary);
};

class TrackingGeometryLambdaVisitor : public TrackingGeometryVisitor {
 public:
  struct Config {
    std::function<void(const TrackingVolume&)> volume{};
    std::function<void(const Portal&)> portal{};
    std::function<void(const Surface&)> surface{};
    std::function<void(const Layer&)> layer{};
    std::function<void(const BoundarySurfaceT<TrackingVolume>&)> boundary{};
  };

  explicit TrackingGeometryLambdaVisitor(Config&& config);

  void visitVolume(const TrackingVolume& volume) override;

  void visitPortal(const Portal& portal) override;

  void visitSurface(const Surface& surface) override;

  void visitLayer(const Layer& layer) override;

  void visitBoundarySurface(
      const BoundarySurfaceT<TrackingVolume>& boundary) override;

 private:
  Config m_config;
};

/// @brief Mutable visitor interface for modifying the tracking geometry hierarchy
///
/// This visitor allows for non-const access to traverse and modify the tracking
/// geometry components. It's used for operations like geometry construction,
/// material decoration, or geometry ID assignment.
class TrackingGeometryMutableVisitor {
 public:
  virtual ~TrackingGeometryMutableVisitor();

  /// @brief Visit and potentially modify a tracking volume
  /// @param volume The tracking volume being visited
  /// @note Called for each volume in the geometry hierarchy during traversal
  virtual void visitVolume(TrackingVolume& volume);

  /// @brief Visit and potentially modify a portal
  /// @param portal The portal being visited
  /// @note Called for each portal encountered during geometry traversal
  virtual void visitPortal(Portal& portal);

  /// @brief Visit and potentially modify a surface
  /// @param surface The surface being visited
  /// @note Called for each surface encountered during geometry traversal
  virtual void visitSurface(Surface& surface);

  // Gen 1
  /// @brief Visit and potentially modify a detector layer
  /// @param layer The layer being visited
  /// @note Called for each layer encountered during geometry traversal
  virtual void visitLayer(Layer& layer);

  /// @brief Visit and potentially modify a boundary surface
  /// @param boundary The boundary surface being visited
  /// @note Called for each boundary surface encountered during geometry traversal
  virtual void visitBoundarySurface(BoundarySurfaceT<TrackingVolume>& boundary);
};

class TrackingGeometryLambdaMutableVisitor
    : public TrackingGeometryMutableVisitor {
 public:
  struct Config {
    std::function<void(TrackingVolume&)> volume{};
    std::function<void(Portal&)> portal{};
    std::function<void(Surface&)> surface{};
    std::function<void(Layer&)> layer{};
    std::function<void(BoundarySurfaceT<TrackingVolume>&)> boundary{};
  };

  explicit TrackingGeometryLambdaMutableVisitor(Config&& config);

  void visitVolume(TrackingVolume& volume) override;

  void visitPortal(Portal& portal) override;

  void visitSurface(Surface& surface) override;

  void visitLayer(Layer& layer) override;

  void visitBoundarySurface(
      BoundarySurfaceT<TrackingVolume>& boundary) override;

 private:
  Config m_config;
};

}  // namespace Acts

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts {

class TrackingVolume;
class Portal;
class Surface;
class Layer;

template <typename T>
class BoundarySurfaceT;

class TrackingGeometryVisitor {
 public:
  virtual ~TrackingGeometryVisitor();

  virtual void visitVolume(const TrackingVolume& volume);

  virtual void visitPortal(const Portal& portla);

  virtual void visitSurface(const Surface& surface);

  // Gen 1
  virtual void visitLayer(const Layer& layer);
  virtual void visitBoundarySurface(
      const BoundarySurfaceT<TrackingVolume>& boundary);
};

class TrackingGeometryMutableVisitor {
 public:
  virtual ~TrackingGeometryMutableVisitor();

  virtual void visitVolume(TrackingVolume&);

  virtual void visitPortal(Portal&);

  virtual void visitSurface(Surface&);

  // Gen 1
  virtual void visitLayer(Layer&);
  virtual void visitBoundarySurface(BoundarySurfaceT<TrackingVolume>& boundary);
};

}  // namespace Acts

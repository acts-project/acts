// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

class IVisualization;

/// Class which implements a cutout cylinder. This shape is bascially a
/// cylinder, with another, smaller cylinder subtracted from the center.
/// --------------------- rmax
/// |                   |
/// |    |---------|    | rmed
/// |    |         |    |
/// ------         ------ rmin
///       -- dz2 --
/// -------- dz1 -------
///
///
class CutoutCylinderVolumeBounds : public VolumeBounds {
 public:
  /// Constructor from defining parameters
  ///
  /// @param rmin Minimum radius at the "choke points"
  /// @param rmed The medium radius (outer radius of the cutout)
  /// @param rmax The outer radius of the overall shape
  /// @param dz1 The longer halflength of the shape
  /// @param dz2 The shorter halflength of the shape
  CutoutCylinderVolumeBounds(double rmin, double rmed, double rmax, double dz1,
                             double dz2)
      : m_rmin(rmin), m_rmed(rmed), m_rmax(rmax), m_dz1(dz1), m_dz2(dz2) {}

  /// Virtual default constructor
  ~CutoutCylinderVolumeBounds() override = default;

  /// Clone method.
  /// @return Pointer to a copy of the shape
  VolumeBounds* clone() const override;

  /// Inside method to test whether a point is inside the shape
  ///
  /// @param gpos The point to test
  /// @param tol The tolerance to test with
  /// @return Whether the point is inside or not.
  bool inside(const Vector3D& gpos, double tol = 0) const override;

  /// Method to decompose the Bounds into Surfaces
  ///
  /// @param transform is the transform to position the surfaces in 3D space
  ///
  /// @return vector of surfaces from the decopmosition
  ///
  std::vector<std::shared_ptr<const Surface>> decomposeToSurfaces(
      const Transform3D* transform = nullptr) const override;

  /// Construct bounding box for this shape
  ///
  /// @param trf Optional transform
  /// @param envelope Optional envelope to add / subtract from min/max
  /// @param entity Entity to associate this bounding box with
  /// @return Constructed bounding box
  Volume::BoundingBox boundingBox(const Transform3D* trf = nullptr,
                                  const Vector3D& envelope = {0, 0, 0},
                                  const Volume* entity = nullptr) const final;

  /// Write information about this instance to an outstream
  ///
  /// @param sl The outstream
  /// @return The outstream
  std::ostream& toStream(std::ostream& sl) const override;

  /// Return the minimum radius
  /// @return The minimum radius
  double rMin() const { return m_rmin; }

  /// Return the medium radius
  /// @return The medium radius
  double rMed() const { return m_rmed; }

  /// Return the maximum radius
  /// @return The maximum radius
  double rMax() const { return m_rmax; }

  /// Return the longer halflength in z.
  /// @return The halflength
  double dZ1() const { return m_dz1; }

  /// Return the shorter halflength in z.
  /// @return The halflength
  double dZ2() const { return m_dz2; }

 private:
  double m_rmin;
  double m_rmed;
  double m_rmax;
  double m_dz1;
  double m_dz2;
};

}  // namespace Acts

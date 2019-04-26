// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <ostream>

#include "Acts/Geometry/Volume.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

class IVisualization;

class GenericCuboidVolumeBounds : public VolumeBounds {
 public:
  GenericCuboidVolumeBounds() = delete;

  /// Constructor from a set of vertices
  /// The ordering is considered to be:
  /// - the first 4 vertices are the "top" face
  /// - the second 4 vertices are the "bottom" face
  /// - both faces are given in counter clock wise order
  GenericCuboidVolumeBounds(const std::array<Acts::Vector3D, 8>& vertices);

  ~GenericCuboidVolumeBounds() override = default;

  ///  clone() method to make deep copy in Volume copy constructor and for
  /// assigment operator  of the Surface class.
  VolumeBounds* clone() const override;

  /// Checking if position given in volume frame is inside
  ///
  /// @param gpos is the global position to be checked
  /// @param tol is the tolerance applied for the inside check
  ///
  /// @return boolean indicating if the position is inside
  bool inside(const Vector3D& gpos, double tol = 0.) const override;

  /// Method to decompose the Bounds into Surfaces
  /// the Volume can turn them into BoundarySurfaces
  ///
  /// @param transform is the 3D transform to be applied to the boundary
  /// surfaces to position them in 3D space
  /// @note this is factory method
  ///
  /// @return a vector of surfaces bounding this volume
  std::vector<std::shared_ptr<const Surface>> decomposeToSurfaces(
      const Transform3D* transform) const override;

  /**
   * Construct bounding box for this shape
   * @param trf Optional transform
   * @param envelope Optional envelope to add / subtract from min/max
   * @param entity Entity to associate this bounding box with
   * @return Constructed bounding box
   */
  Volume::BoundingBox boundingBox(const Transform3D* trf = nullptr,
                                  const Vector3D& envelope = {0, 0, 0},
                                  const Volume* entity = nullptr) const final;

  ///
  /// @param sl is the output stream to be written into
  std::ostream& toStream(std::ostream& sl) const override;

  /**
   * Draw this shape using a visualization helper
   * @param helper The visualizatin helper
   * @param transform Optional transformation matrix
   */
  void draw(IVisualization& helper,
            const Transform3D& transform = Transform3D::Identity()) const;

 private:
  std::array<Vector3D, 8> m_vertices;
  std::array<Vector3D, 6> m_normals;
};
}  // namespace Acts

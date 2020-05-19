// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Volume.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"

#include <array>
#include <ostream>
#include <vector>

namespace Acts {

class IVisualization;

class GenericCuboidVolumeBounds : public VolumeBounds {
 public:
  static constexpr size_t eSize = 24;

  GenericCuboidVolumeBounds() = delete;

  /// Constructor from a set of vertices
  ///
  /// @param vertices The set of input vertices
  ///
  /// The ordering is considered to be:
  /// - the first 4 vertices are the "top" face
  /// - the second 4 vertices are the "bottom" face
  /// - both faces are given in counter clock wise order
  GenericCuboidVolumeBounds(
      const std::array<Acts::Vector3D, 8>& vertices) noexcept(false);

  /// Constructor from a fixed size array
  ///
  /// @param values The input values
  GenericCuboidVolumeBounds(const std::array<double, eSize>& values) noexcept(
      false);

  ~GenericCuboidVolumeBounds() override = default;

  VolumeBounds::BoundsType type() const final {
    return VolumeBounds::eGenericCuboid;
  }

  /// Return the bound values as dynamically sized vector
  ///
  /// @return this returns a copy of the internal values
  std::vector<double> values() const final;

  /// Checking if position given in volume frame is inside
  ///
  /// @param gpos is the global position to be checked
  /// @param tol is the tolerance applied for the inside check
  ///
  /// @return boolean indicating if the position is inside
  bool inside(const Vector3D& gpos, double tol = 0.) const override;

  /// Oriented surfaces, i.e. the decomposed boundary surfaces and the
  /// according navigation direction into the volume given the normal
  /// vector on the surface
  ///
  /// @param transform is the 3D transform to be applied to the boundary
  /// surfaces to position them in 3D space
  ///
  /// It will throw an exception if the orientation prescription is not adequate
  ///
  /// @return a vector of surfaces bounding this volume
  OrientedSurfaces orientedSurfaces(
      const Transform3D* transform = nullptr) const override;

  /// Construct bounding box for this shape
  /// @param trf Optional transform
  /// @param envelope Optional envelope to add / subtract from min/max
  /// @param entity Entity to associate this bounding box with
  /// @return Constructed bounding box
  Volume::BoundingBox boundingBox(const Transform3D* trf = nullptr,
                                  const Vector3D& envelope = {0, 0, 0},
                                  const Volume* entity = nullptr) const final;

  /// @param sl is the output stream to be written into
  std::ostream& toStream(std::ostream& sl) const override;

  /// Draw this shape using a visualization helper
  /// @param helper The visualizatin helper
  /// @param transform Optional transformation matrix
  ///
  void draw(IVisualization& helper,
            const Transform3D& transform = Transform3D::Identity()) const;

 private:
  std::array<Vector3D, 8> m_vertices;
  std::array<Vector3D, 6> m_normals;

  /// Private helper method to contruct the Volume bounds
  /// to be called by the constructors, from the ordered input vertices
  void construct() noexcept(false);
};

}  // namespace Acts

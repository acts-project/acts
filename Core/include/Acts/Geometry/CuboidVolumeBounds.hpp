// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/BoundingBox.hpp"

#include <array>
#include <cmath>
#include <iomanip>
#include <iosfwd>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace Acts {

class RectangleBounds;

/// @class CuboidVolumeBounds
///
/// Bounds for a cubical Volume, the orientedSurfaces(...) method creates a
/// vector of 6 surfaces:
///
///  BoundarySurfaceFace [index]:
///
///    - negativeFaceXY [0] : Rectangular Acts::PlaneSurface, parallel to \f$ xy
/// \f$ plane at negative \f$ z \f$
///    - positiveFaceXY [1] : Rectangular Acts::PlaneSurface, parallel to \f$ xy
/// \f$ plane at positive \f$ z \f$
///    - negativeFaceXY [2] : Rectangular Acts::PlaneSurface, attached to \f$ yz
/// \f$ plane at negative \f$ x \f$
///    - positiveFaceXY [3] : Rectangular Acts::PlaneSurface, attached to \f$ yz
/// \f$ plane at negative \f$ x \f$
///    - negativeFaceXY [4] : Rectangular Acts::PlaneSurface, parallel to \f$ zx
/// \f$ plane at negative \f$ y \f$
///    - positiveFaceXY [5] : Rectangular Acts::PlaneSurface, parallel to \f$ zx
/// \f$ plane at positive \f$ y \f$
///
class CuboidVolumeBounds : public VolumeBounds {
 public:
  /// @enum BoundValues for streaming and access
  enum BoundValues : unsigned int {
    eHalfLengthX = 0,
    eHalfLengthY = 1,
    eHalfLengthZ = 2,
    eSize
  };

  CuboidVolumeBounds() = delete;

  /// Constructor - the box boundaries
  ///
  /// @param halex is the half length of the cube in x
  /// @param haley is the half length of the cube in y
  /// @param halez is the half length of the cube in z
  CuboidVolumeBounds(double halex, double haley, double halez) noexcept(false);

  /// Constructor - from a fixed size array
  ///
  /// @param values iw the bound values
  CuboidVolumeBounds(const std::array<double, eSize>& values) noexcept(false)
      : m_values(values) {
    checkConsistency();
    buildSurfaceBounds();
  }

  /// Copy Constructor
  ///
  /// @param bobo is the source volume bounds to be copied
  CuboidVolumeBounds(const CuboidVolumeBounds& bobo);

  /// Assignment operator
  ///
  /// @param bobo is the source volume bounds to be assigned
  CuboidVolumeBounds& operator=(const CuboidVolumeBounds& bobo);

  ~CuboidVolumeBounds() override = default;

  VolumeBounds::BoundsType type() const final { return VolumeBounds::eCuboid; }

  /// Return the bound values as dynamically sized vector
  ///
  /// @return this returns a copy of the internal values
  std::vector<double> values() const final;

  /// This method checks if position in the 3D volume
  /// frame is inside the cylinder
  ///
  /// @param pos is the position in volume frame to be checked
  /// @param tol is the absolute tolerance to be applied
  bool inside(const Vector3& pos, double tol = 0.) const override;

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
      const Transform3& transform = Transform3::Identity()) const override;

  /// Construct bounding box for this shape
  /// @param trf Optional transform
  /// @param envelope Optional envelope to add / subtract from min/max
  /// @param entity Entity to associate this bounding box with
  /// @return Constructed bounding box
  Volume::BoundingBox boundingBox(const Transform3* trf = nullptr,
                                  const Vector3& envelope = {0, 0, 0},
                                  const Volume* entity = nullptr) const final;

  /// Binning borders in double
  ///
  /// @param bValue is the binning schema used
  ///
  /// @return float offset to be used for the binning
  double binningBorder(BinningValue bValue) const final;

  /// Access to the bound values
  /// @param bValue the class nested enum for the array access
  double get(BoundValues bValue) const { return m_values[bValue]; }

  /// Output Method for std::ostream
  ///
  /// @param sl is ostream operator to be dumped into
  std::ostream& toStream(std::ostream& sl) const override;

 private:
  /// Templated dumpT method
  /// @tparam stream_t The type for the dump stream
  /// @param dt The dump stream object
  template <class stream_t>
  stream_t& dumpT(stream_t& dt) const;

  /// The bound values ordered in a fixed size array
  std::array<double, eSize> m_values;

  std::shared_ptr<const RectangleBounds> m_xyBounds{nullptr};
  std::shared_ptr<const RectangleBounds> m_yzBounds{nullptr};
  std::shared_ptr<const RectangleBounds> m_zxBounds{nullptr};

  /// Create the surface bounds
  void buildSurfaceBounds();

  /// Check the input values for consistency,
  /// will throw a logic_exception if consistency is not given
  void checkConsistency() noexcept(false);
};

inline bool CuboidVolumeBounds::inside(const Vector3& pos, double tol) const {
  return (std::abs(pos.x()) <= get(eHalfLengthX) + tol &&
          std::abs(pos.y()) <= get(eHalfLengthY) + tol &&
          std::abs(pos.z()) <= get(eHalfLengthZ) + tol);
}

inline std::vector<double> CuboidVolumeBounds::values() const {
  std::vector<double> valvector;
  valvector.insert(valvector.begin(), m_values.begin(), m_values.end());
  return valvector;
}

inline void CuboidVolumeBounds::checkConsistency() noexcept(false) {
  if (get(eHalfLengthX) <= 0 || get(eHalfLengthY) <= 0 ||
      get(eHalfLengthZ) <= 0.) {
    throw std::invalid_argument(
        "CuboidVolumeBounds: invalid input, zero or negative.");
  }
}

template <class stream_t>
stream_t& CuboidVolumeBounds::dumpT(stream_t& dt) const {
  dt << std::setiosflags(std::ios::fixed);
  dt << std::setprecision(5);
  dt << "Acts::CuboidVolumeBounds: (halfLengthX, halfLengthY, halfLengthZ) = ";
  dt << "(" << get(eHalfLengthX) << ", " << get(eHalfLengthY) << ", "
     << get(eHalfLengthZ) << ")";
  return dt;
}
}  // namespace Acts

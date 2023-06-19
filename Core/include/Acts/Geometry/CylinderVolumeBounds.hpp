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
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <array>
#include <cmath>
#include <iomanip>
#include <iosfwd>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace Acts {

class CylinderBounds;
class RadialBounds;
class PlanarBounds;

/// @class CylinderVolumeBounds
///
/// Bounds for a cylindrical Volume, the orientedSurfaces(..) method creates a
/// vector of up to 6 surfaces:
///
/// case A) 3 Surfaces (full cylindrical tube):
///  BoundarySurfaceFace [index]:
///  - negativeFaceXY [0] : Acts::DiscSurface with \f$ r_{inner}=0 \f$,
///                         parallel to \f$ xy \f$ plane at negative \f$ z\f$
///  - positiveFaceXY [1] : Acts::DiscSurface with \f$ r_{inner}=0 \f$,
///                         parallel to \f$ xy \f$ plane at positive \f$ z\f$
///  - cylinderCover  [2] : Acts::CylinderSurface confining the Acts::Volume
///
/// case B) 4 Surfaces (tube with inner and outer radius):
///  BoundarySurfaceFace [index]:
///  - negativeFaceXY [0] : Acts::DiscSurface with \f$ r_{inner}>0 \f$,
///                         parallel to \f$ xy \f$ plane at negative \f$ z\f$
///  - positiveFaceXY [1] : Acts::DiscSurface with \f$ r_{inner}>0 \f$,
///                         parallel to \f$ xy \f$ plane at positive \f$ z\f$
///  - tubeOuterCover [2] : Acts::CylinderSurface with \f$ r = r_{outer} \f$
///  - tubeInnerCover [3] : Acts::CylinderSurface with \f$ r = r_{inner} \f$
///
/// case C) 6 Surfaces (sectoral tube with inner and outer radius):
///  BoundarySurfaceFace [index]:
///   - negativeFaceXY  [0] : Acts::DiscSurface with \f$ r_{inner}>0\f$
///                           and \f$ \phi < \pi \f$,
///                           parallel to \f$ xy \f$ plane at negative \f$z\f$
///   - positiveFaceXY  [1] : Acts::DiscSurface with \f$ r_{inner}>0 \f$
///                           and \f$ \phi < \pi \f$,
///                           parallel to \f$ xy \f$ plane at positive \f$z\f$
///   - tubeSectorOuterCover  [2] : Acts::CylinderSurface with
///                                 \f$ r = r_{outer}\f$
///   - tubeSectorInnerCover  [3] : Acts::CylinderSurface with
///                                 \f$ r = r_{inner} \f$
///   - tubeSectorNegativePhi [4] : Rectangular Acts::PlaneSurface attached to
///                 [0] and [1] at negative \f$ \phi \f$
///                      - tubeSectorNegativePhi [5] :
//                          Rectangular Acts::PlaneSurface attached to
///                 [0] and [1] at positive \f$ \phi \f$
///
class CylinderVolumeBounds : public VolumeBounds {
 public:
  /// @enum BoundValues for streaming and access
  enum BoundValues : unsigned int {
    eMinR = 0,
    eMaxR = 1,
    eHalfLengthZ = 2,
    eHalfPhiSector = 3,
    eAveragePhi = 4,
    eBevelMinZ = 5,
    eBevelMaxZ = 6,
    eSize
  };

  CylinderVolumeBounds() = delete;

  /// Constructor
  ///
  /// @param rmin The inner radius of the cylinder
  /// @param rmax The outer radius of the cylinder
  /// @param halfz The half length in z
  /// @param halfphi The half lopening angle
  /// @param avgphi The average phi value
  /// @param bevelMinZ The bevel angle, in radians, for the negative side
  /// @param bevelMaxZ The bevel angle, in radians, for the positive side
  CylinderVolumeBounds(double rmin, double rmax, double halfz,
                       double halfphi = M_PI, double avgphi = 0.,
                       double bevelMinZ = 0.,
                       double bevelMaxZ = 0.) noexcept(false)
      : m_values({rmin, rmax, halfz, halfphi, avgphi, bevelMinZ, bevelMaxZ}) {
    checkConsistency();
    buildSurfaceBounds();
  }

  /// Constructor - from a fixed size array
  ///
  /// @param values The bound values
  CylinderVolumeBounds(const std::array<double, eSize>& values) noexcept(false)
      : m_values(values) {
    checkConsistency();
    buildSurfaceBounds();
  }

  /// Constructor - extruded from cylinder bounds and thickness
  ///
  /// @param cBounds the cylinder bounds
  /// @param thickness of the extrusion
  CylinderVolumeBounds(const CylinderBounds& cBounds,
                       double thickness) noexcept(false);

  /// Constructor - extruded from radial bounds and thickness
  ///
  /// @param rBounds the Radial bounds
  /// @param thickness
  CylinderVolumeBounds(const RadialBounds& rBounds,
                       double thickness) noexcept(false);

  /// Copy Constructor
  ///
  /// @param cylbo is the source cylinder volume bounds for the copy
  CylinderVolumeBounds(const CylinderVolumeBounds& cylbo) = default;

  ~CylinderVolumeBounds() override = default;
  CylinderVolumeBounds& operator=(const CylinderVolumeBounds& cylbo) = default;

  VolumeBounds::BoundsType type() const final {
    return VolumeBounds::eCylinder;
  }

  /// Return the bound values as dynamically sized vector
  ///
  /// @return this returns a copy of the internal values
  std::vector<double> values() const final;

  /// This method checks if position in the 3D volume
  /// frame is inside the cylinder
  ///
  /// @param pos is a global position to be checked
  /// @param tol is the tolerance for the check
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

  /// Binning offset - overloaded for some R-binning types
  ///
  /// @param bValue is the type used for the binning
  Vector3 binningOffset(BinningValue bValue) const override;

  /// Binning borders in double
  ///
  /// @param bValue is the type used for the binning
  double binningBorder(BinningValue bValue) const override;

  /// Output Method for std::ostream
  std::ostream& toStream(std::ostream& sl) const override;

  /// Access to the bound values
  /// @param bValue the class nested enum for the array access
  double get(BoundValues bValue) const { return m_values[bValue]; }

 private:
  /// The internal version of the bounds can be float/double
  std::array<double, eSize> m_values{};
  /// Bounds of the inner CylinderBounds
  std::shared_ptr<const CylinderBounds> m_innerCylinderBounds{nullptr};
  /// Bounds of the inner CylinderBounds
  std::shared_ptr<const CylinderBounds> m_outerCylinderBounds{nullptr};
  /// Bounds of the bottom/top Radial
  std::shared_ptr<const RadialBounds> m_discBounds{nullptr};
  /// Bounds of the sector planes
  std::shared_ptr<const PlanarBounds> m_sectorPlaneBounds{nullptr};

  /// Check the input values for consistency,
  /// will throw a logic_exception if consistency is not given
  void checkConsistency() noexcept(false);

  /// Helper method to create the surface bounds
  void buildSurfaceBounds();

  /// Templated dumpT method
  /// @tparam stream_t The type fo the dump stream
  /// @param dt The dump stream object
  template <class stream_t>
  stream_t& dumpT(stream_t& dt) const;
};

inline bool CylinderVolumeBounds::inside(const Vector3& pos, double tol) const {
  using VectorHelpers::perp;
  using VectorHelpers::phi;
  double ros = perp(pos);
  bool insidePhi = cos(phi(pos)) >= cos(get(eHalfPhiSector)) - tol;
  bool insideR = insidePhi
                     ? ((ros >= get(eMinR) - tol) && (ros <= get(eMaxR) + tol))
                     : false;
  bool insideZ =
      insideR ? (std::abs(pos.z()) <= get(eHalfLengthZ) + tol) : false;
  return (insideZ && insideR && insidePhi);
}

inline Vector3 CylinderVolumeBounds::binningOffset(BinningValue bValue)
    const {  // the medium radius is taken for r-type binning
  if (bValue == Acts::binR || bValue == Acts::binRPhi) {
    return Vector3(0.5 * (get(eMinR) + get(eMaxR)), 0., 0.);
  }
  return VolumeBounds::binningOffset(bValue);
}

inline double CylinderVolumeBounds::binningBorder(BinningValue bValue) const {
  if (bValue == Acts::binR) {
    return 0.5 * (get(eMaxR) - get(eMinR));
  }
  if (bValue == Acts::binZ) {
    return get(eHalfLengthZ);
  }
  return VolumeBounds::binningBorder(bValue);
}

template <class stream_t>
stream_t& CylinderVolumeBounds::dumpT(stream_t& dt) const {
  dt << std::setiosflags(std::ios::fixed);
  dt << std::setprecision(5);
  dt << "Acts::CylinderVolumeBounds: (rMin, rMax, halfZ, halfPhi, "
        "averagePhi, minBevelZ, maxBevelZ) = ";
  dt << get(eMinR) << ", " << get(eMaxR) << ", " << get(eHalfLengthZ) << ", "
     << get(eHalfPhiSector) << ", " << get(eAveragePhi) << ", "
     << get(eBevelMinZ) << ", " << get(eBevelMaxZ);
  return dt;
}

inline std::vector<double> CylinderVolumeBounds::values() const {
  std::vector<double> valvector;
  valvector.insert(valvector.begin(), m_values.begin(), m_values.end());
  return valvector;
}

inline void CylinderVolumeBounds::checkConsistency() noexcept(false) {
  if (get(eMinR) < 0. or get(eMaxR) <= 0. or get(eMinR) >= get(eMaxR)) {
    throw std::invalid_argument("CylinderVolumeBounds: invalid radial input.");
  }
  if (get(eHalfLengthZ) <= 0) {
    throw std::invalid_argument(
        "CylinderVolumeBounds: invalid longitudinal input.");
  }
  if (get(eHalfPhiSector) < 0. or get(eHalfPhiSector) > M_PI) {
    throw std::invalid_argument(
        "CylinderVolumeBounds: invalid phi sector setup.");
  }
  if (get(eAveragePhi) != detail::radian_sym(get(eAveragePhi))) {
    throw std::invalid_argument(
        "CylinderVolumeBounds: invalid phi positioning.");
  }
  if (get(eBevelMinZ) != detail::radian_sym(get(eBevelMinZ))) {
    throw std::invalid_argument("CylinderBounds: invalid bevel at min Z.");
  }
  if (get(eBevelMaxZ) != detail::radian_sym(get(eBevelMaxZ))) {
    throw std::invalid_argument("CylinderBounds: invalid bevel at max Z.");
  }
}

}  // namespace Acts

// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <array>
#include <cmath>
#include <exception>
#include <vector>

namespace Acts {

class Surface;
class CylinderBounds;
class RadialBounds;
class PlanarBounds;
class IVisualization;

/// @class CylinderVolumeBounds
///
/// Bounds for a cylindrical Volume, the decomposeToSurfaces method creates a
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
///  @image html CylinderVolumeBounds_decomp.gif

class CylinderVolumeBounds : public VolumeBounds {
 public:
  /// @enum BoundValues for streaming and access
  enum BoundValues {
    eMinR = 0,
    eMaxR = 1,
    eHalfLengthZ = 2,
    eHalfPhiSector = 3,
    eAveragePhi = 4,
    eSize = 5
  };

  CylinderVolumeBounds() = delete;

  /// Constructor
  ///
  /// @param rmin The inner radius of the cylinder
  /// @param rmax The outer radius of the cylinder
  /// @param halfz The half length in z
  /// @param halfphi The half lopening angle
  /// @param avgphi The average phi value
  CylinderVolumeBounds(double rmin, double rmax, double halfz,
                       double halfphi = M_PI,
                       double avgphi = 0.) noexcept(false)
      : m_values({rmin, rmax, halfz, halfphi, avgphi}) {
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
  /// @param cbounds the cylinder bounds
  /// @param thickness of the extrusion
  CylinderVolumeBounds(const CylinderBounds& cBounds,
                       double thickness) noexcept(false);

  /// Constructor - extruded from radial bounds and thickness
  ///
  /// @param rbounds the Radial bounds
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
  bool inside(const Vector3D& pos, double tol = 0.) const override;

  /// Method to decompose the Bounds into boundarySurfaces
  /// @param transformPtr is the transform where the boundary surfaces are
  /// situated
  std::vector<std::shared_ptr<const Surface>> decomposeToSurfaces(
      const Transform3D* transformPtr = nullptr) const override;

  /// Construct bounding box for this shape
  /// @param trf Optional transform
  /// @param envelope Optional envelope to add / subtract from min/max
  /// @param entity Entity to associate this bounding box with
  /// @return Constructed bounding box
  Volume::BoundingBox boundingBox(const Transform3D* trf = nullptr,
                                  const Vector3D& envelope = {0, 0, 0},
                                  const Volume* entity = nullptr) const final;

  /// Binning offset - overloaded for some R-binning types
  ///
  /// @param bValue is the type used for the binning
  Vector3D binningOffset(BinningValue bValue) const override;

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
  std::array<double, eSize> m_values;
  /// Bounds of the inner CylinderSurfaces
  std::shared_ptr<const CylinderBounds> m_innerCylinderBounds{nullptr};
  /// Bounds of the inner CylinderSurfaces
  std::shared_ptr<const CylinderBounds> m_outerCylinderBounds{nullptr};
  /// Bounds of the bottom/top DiscSurface
  std::shared_ptr<const RadialBounds> m_discBounds{nullptr};
  /// Bounds of the sector planes
  std::shared_ptr<const PlanarBounds> m_sectorPlaneBounds{nullptr};

  /// Check the input values for consistency,
  /// will throw a logic_exception if consistency is not given
  void checkConsistency() noexcept(false);

  /// Helper method to create the surface bounds
  void buildSurfaceBounds();

  /// templated dumpT method
  template <class T>
  T& dumpT(T& tstream) const;
};

inline bool CylinderVolumeBounds::inside(const Vector3D& pos,
                                         double tol) const {
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

inline Vector3D CylinderVolumeBounds::binningOffset(BinningValue bValue)
    const {  // the medium radius is taken for r-type binning
  if (bValue == Acts::binR || bValue == Acts::binRPhi) {
    return Vector3D(0.5 * (get(eMinR) + get(eMaxR)), 0., 0.);
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

template <class T>
T& CylinderVolumeBounds::dumpT(T& tstream) const {
  tstream << std::setiosflags(std::ios::fixed);
  tstream << std::setprecision(5);
  tstream << "Acts::CylinderVolumeBounds: (rMin, rMax, halfZ, halfPhi, "
             "averagePhi) = ";
  tstream << get(eMinR) << ", " << get(eMaxR) << ", " << get(eHalfLengthZ)
          << ", " << get(eHalfPhiSector) << get(eAveragePhi);
  return tstream;
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
}

}  // namespace Acts
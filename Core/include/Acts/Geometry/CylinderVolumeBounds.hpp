// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <cmath>
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {

class Surface;
class CylinderBounds;
class DiscBounds;
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
  /// @enum BoundValues for readability
  enum BoundValues {
    bv_innerRadius = 0,
    bv_outerRadius = 1,
    bv_halfPhiSector = 2,
    bv_halfZ = 3,
    bv_length = 4
  };

  /// Default Constructor
  CylinderVolumeBounds();

  /// Constructor - full cylinder
  ///
  /// @param radius is the outer radius of the cylinder
  /// @param halez is the half length in z
  CylinderVolumeBounds(double radius, double halez);

  /// Constructor - extruded cylinder
  ///
  /// @param rinner is the inner radius of the cylinder
  /// @param router is the outer radius of the cylinder
  /// @param halez is the half length in z
  CylinderVolumeBounds(double rinner, double router, double halez);

  /// Constructor - extruded cylinder
  ///
  /// @param rinner is the inner radius of the cylinder
  /// @param router is the outer radius of the cylinder
  /// @param haphi is the half opening angle
  /// @param halez is the half length in z
  CylinderVolumeBounds(double rinner, double router, double haphi,
                       double halez);

  /// Constructor - from cylinder bounds and thickness
  ///
  /// @param cbounds the cylinder bounds
  /// @param thickness
  CylinderVolumeBounds(const CylinderBounds& cBounds, double thickness);

  /// Constructor - from radial bounds and thickness
  ///
  /// @param rbounds the Radial bounds
  /// @param thickness
  CylinderVolumeBounds(const RadialBounds& rBounds, double thickness);

  /// Copy Constructor
  ///
  /// @param cylbo is the source cylinder volume bounds for the copy
  CylinderVolumeBounds(const CylinderVolumeBounds& cylbo);

  /// Destructor
  ~CylinderVolumeBounds() override;

  /// Assignment operator
  CylinderVolumeBounds& operator=(const CylinderVolumeBounds& cylbo);

  /// Virtual constructor
  CylinderVolumeBounds* clone() const override;

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

  /// This method returns the inner radius
  double innerRadius() const;

  /// This method returns the outer radius
  double outerRadius() const;

  /// This method returns the medium radius
  double mediumRadius() const;

  /// This method returns the delta radius
  double deltaRadius() const;

  /// This method returns the halfPhiSector angle
  double halfPhiSector() const;

  /// This method returns the halflengthZ
  double halflengthZ() const;

  /// Output Method for std::ostream
  std::ostream& toStream(std::ostream& sl) const override;

 private:
  /// templated dumpT method
  template <class T>
  T& dumpT(T& tstream) const;

  /// This method returns the associated CylinderBounds
  /// of the inner CylinderSurfaces.
  std::shared_ptr<const CylinderBounds> innerCylinderBounds() const;

  /// This method returns the associated CylinderBounds
  /// of the inner CylinderSurfaces.
  std::shared_ptr<const CylinderBounds> outerCylinderBounds() const;

  /// This method returns the associated RadialBounds
  /// for the bottom/top DiscSurface
  std::shared_ptr<const DiscBounds> discBounds() const;

  /// This method returns the associated PlaneBounds
  /// limiting a sectoral CylinderVolume
  std::shared_ptr<const PlanarBounds> sectorPlaneBounds() const;

  /// The internal version of the bounds can be float/double
  std::vector<TDD_real_t> m_valueStore;

  /// numerical stability
  /// @todo unify the numerical stability checks
  static const double s_numericalStable;
};

inline CylinderVolumeBounds* CylinderVolumeBounds::clone() const {
  return new CylinderVolumeBounds(*this);
}

inline bool CylinderVolumeBounds::inside(const Vector3D& pos,
                                         double tol) const {
  using VectorHelpers::perp;
  using VectorHelpers::phi;
  double ros = perp(pos);
  bool insidePhi = cos(phi(pos)) >= cos(m_valueStore[bv_halfPhiSector]) - tol;
  bool insideR = insidePhi ? ((ros >= m_valueStore[bv_innerRadius] - tol) &&
                              (ros <= m_valueStore[bv_outerRadius] + tol))
                           : false;
  bool insideZ =
      insideR ? (std::abs(pos.z()) <= m_valueStore[bv_halfZ] + tol) : false;
  return (insideZ && insideR && insidePhi);
}

inline Vector3D CylinderVolumeBounds::binningOffset(BinningValue bValue)
    const {  // the medium radius is taken for r-type binning
  if (bValue == Acts::binR || bValue == Acts::binRPhi) {
    return Vector3D(mediumRadius(), 0., 0.);
  }
  return VolumeBounds::binningOffset(bValue);
}

inline double CylinderVolumeBounds::binningBorder(BinningValue bValue)
    const {  // the medium radius is taken for r-type binning
  if (bValue == Acts::binR) {
    return 0.5 * deltaRadius();
  }
  if (bValue == Acts::binZ) {
    return halflengthZ();
  }
  return VolumeBounds::binningBorder(bValue);
}

inline double CylinderVolumeBounds::innerRadius() const {
  return m_valueStore.at(bv_innerRadius);
}

inline double CylinderVolumeBounds::outerRadius() const {
  return m_valueStore.at(bv_outerRadius);
}

inline double CylinderVolumeBounds::mediumRadius() const {
  return 0.5 *
         (m_valueStore.at(bv_innerRadius) + m_valueStore.at(bv_outerRadius));
}

inline double CylinderVolumeBounds::deltaRadius() const {
  return (m_valueStore.at(bv_outerRadius) - m_valueStore.at(bv_innerRadius));
}

inline double CylinderVolumeBounds::halfPhiSector() const {
  return m_valueStore.at(bv_halfPhiSector);
}

inline double CylinderVolumeBounds::halflengthZ() const {
  return m_valueStore.at(bv_halfZ);
}

template <class T>
T& CylinderVolumeBounds::dumpT(T& tstream) const {
  tstream << std::setiosflags(std::ios::fixed);
  tstream << std::setprecision(5);
  tstream << "Acts::CylinderVolumeBounds: (rMin, rMax, halfPhi, halfZ) = ";
  tstream << m_valueStore.at(bv_innerRadius) << ", "
          << m_valueStore.at(bv_outerRadius) << ", "
          << m_valueStore.at(bv_halfPhiSector) << ", "
          << m_valueStore.at(bv_halfZ);
  return tstream;
}
}  // namespace Acts
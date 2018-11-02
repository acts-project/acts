// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrapezoidVolumeBounds.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Volumes/VolumeBounds.hpp"

namespace Acts {

class Surface;
class RectangleBounds;
class TrapezoidBounds;

/// @class TrapezoidVolumeBounds
///
/// Bounds for a trapezoidal shaped Volume, the decomposeToSurfaces method
/// creates a vector of 6 surfaces:
///
///  BoundarySurfaceFace [index]:
///
///  - negativeFaceXY     [0] : Trazpezoidal Acts::PlaneSurface,
///                             parallel to \f$ xy \f$ plane at negative \f$z\f$
///  - positiveFaceXY     [1] : Trazpezoidal Acts::PlaneSurface,
///                             parallel to \f$ xy \f$ plane at positive \f$z\f$
///  - trapezoidFaceAlpha [2] : Rectangular  Acts::PlaneSurface,
///                             attached to [0] and [1] at negative \f$x\f$
///  (associated to alpha)
///  - trapezoidFaceBeta  [3] : Rectangular  Acts::PlaneSurface,
///                             attached to [0] and [1] at positive \f$ x\f$
///  (associated to beta)
///  - negativeFaceZX     [4] : Rectangular  Acts::PlaneSurface,
///                             parallel to \f$ zx \f$ plane at negative \f$y\f$
///  - positiveFaceZX     [5] : Rectangular  Acts::PlaneSurface,
///                             parallel to \f$ zx \f$ plane at positive \f$y\f$
///
///  @image html TrapezoidVolumeBounds_decomp.gif

class TrapezoidVolumeBounds : public VolumeBounds
{
public:
  /// @enum BoundValues for readability
  enum BoundValues {
    bv_minHalfX = 0,  //!< minimal halflength in x
    bv_maxHalfX = 1,  //!< maximal halflength in x
    bv_halfY    = 2,  //!< halflength in y
    bv_halfZ    = 3,  //!< halflength in z
    bv_alpha    = 4,  //!< opening angle alpha (in point A)
    bv_beta     = 5,  //!< opening angle beta  (in point B)
    bv_length   = 6   // length of the bounds vector

  };

  /// Default Constructor
  TrapezoidVolumeBounds();

  /// Constructor - the trapezoid boundaries (symmetric trapezoid)
  ///
  /// @param minhalex is the half length in x at minimal y
  /// @param maxhalex is the half length in x at maximal y
  /// @param haley is the half length in y
  /// @param halez is the half length in z
  TrapezoidVolumeBounds(double minhalex,
                        double maxhalex,
                        double haley,
                        double halez);

  /// Constructor - the trapezoid boundaries (arbitrary trapezoid)
  ///
  /// @param minhalex is the half length in x at minimal y
  /// @param haley is the half length in y
  /// @param halez is the half length in z
  /// @param alpha is the openeing angle at -x,-y
  /// @param beta is the openeing angle at +x,-y
  TrapezoidVolumeBounds(double minhalex,
                        double haley,
                        double halez,
                        double alpha,
                        double beta);

  /// Copy Constructor
  /// @param trabo The object to be copied
  TrapezoidVolumeBounds(const TrapezoidVolumeBounds& trabo);

  /// Destructor
  ~TrapezoidVolumeBounds() override;

  /// Assignment operator
  /// @param trabo The object to be assigned
  TrapezoidVolumeBounds&
  operator=(const TrapezoidVolumeBounds& trabo);

  /// Virtual constructor
  TrapezoidVolumeBounds*
  clone() const override;

  /// This method checks if position in the 3D volume frame
  /// is inside the cylinder
  ///
  /// @param pos is the global position to be checked
  /// @param tol is the tolerance applied
  ///
  /// @return boolean indicator if position is inside
  bool
  inside(const Vector3D& pos, double tol = 0.) const override;

  /// Method to decompose the Bounds into Surfaces
  ///
  /// @param transformPtr is the transform to position the surfaces in 3D space
  /// @note this is a factory method
  ///
  /// @return vector of surfaces from the decopmosition
  std::vector<std::shared_ptr<const Surface>>
  decomposeToSurfaces(
      std::shared_ptr<const Transform3D> transformPtr) const override;

  /// This method returns the minimal halflength in local x
  double
  minHalflengthX() const;

  /// This method returns the maximal halflength in local x
  double
  maxHalflengthX() const;

  /// This method returns the halflength in local y
  double
  halflengthY() const;

  /// This method returns the halflength in local z
  double
  halflengthZ() const;

  /// This method returns the opening angle in point A (negative local x)
  double
  alpha() const;

  /// This method returns the opening angle in point B (negative local x)
  double
  beta() const;

  /// Output Method for std::ostream
  std::ostream&
  dump(std::ostream& sl) const override;

private:
  /// Templated dump methos
  template <class T>
  T&
  dumpT(T& dt) const;

  /// This method returns the associated TrapezoidBounds
  /// of the face PlaneSurface parallel to local xy plane
  TrapezoidBounds*
  faceXYTrapezoidBounds() const;

  /// This method returns the associated RecangleBounds
  /// of the face PlaneSurface attached to alpha (negative local x)
  RectangleBounds*
  faceAlphaRectangleBounds() const;

  /// This method returns the associated RectangleBounds
  /// of the face PlaneSurface attached to beta (positive local x)
  RectangleBounds*
  faceBetaRectangleBounds() const;

  /// This method returns the associated RectangleBounds
  /// of the face PlaneSurface parallel to local zx plane, negative local y
  RectangleBounds*
  faceZXRectangleBoundsBottom() const;

  /// This method returns the associated RectangleBounds
  /// of the face PlaneSurface parallel to local zx plane, positive local y
  RectangleBounds*
  faceZXRectangleBoundsTop() const;

  /// the bounds values
  std::vector<TDD_real_t> m_valueStore;
};

inline TrapezoidVolumeBounds*
TrapezoidVolumeBounds::clone() const
{
  return new TrapezoidVolumeBounds(*this);
}

inline double
TrapezoidVolumeBounds::minHalflengthX() const
{
  return m_valueStore.at(bv_minHalfX);
}
inline double
TrapezoidVolumeBounds::maxHalflengthX() const
{
  return m_valueStore.at(bv_maxHalfX);
}
inline double
TrapezoidVolumeBounds::halflengthY() const
{
  return m_valueStore.at(bv_halfY);
}
inline double
TrapezoidVolumeBounds::halflengthZ() const
{
  return m_valueStore.at(bv_halfZ);
}
inline double
TrapezoidVolumeBounds::alpha() const
{
  return m_valueStore.at(bv_alpha);
}
inline double
TrapezoidVolumeBounds::beta() const
{
  return m_valueStore.at(bv_beta);
}

template <class T>
T&
TrapezoidVolumeBounds::dumpT(T& dt) const
{
  dt << std::setiosflags(std::ios::fixed);
  dt << std::setprecision(5);
  dt << "Acts::TrapezoidVolumeBounds: (minhalfX, halfY, halfZ, alpha, beta) = ";
  dt << "(" << m_valueStore.at(bv_minHalfX) << ", " << m_valueStore.at(bv_halfY)
     << ", " << m_valueStore.at(bv_halfZ);
  dt << ", " << m_valueStore.at(bv_alpha) << ", " << m_valueStore.at(bv_beta)
     << ")";
  return dt;
}
}
// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BoundaryCheck.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_BOUNDARYCHECK_H
#define ACTS_SURFACES_BOUNDARYCHECK_H 1

#include <cmath>
#include <limits>
#include <vector>
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/ParameterDefinitions.hpp"

namespace Acts {

/// @struct KDOP
/// maint struct for comparing the extent of arbitrary convex polygons relative
/// to prefixed axes
struct KDOP
{
  float min;
  // Minimum distance (from origin) along axes
  float max;
  // Maximum distance (from origin) along axes
};

///  @class BoundaryCheck
///
///  The BoundaryCheck class allows to steer the way surface
///  boundaries are used for inside/outside checks of
///  parameters.
///
///  These checks are performed in the local 2D frame of the
///  surface and can either be:
///  - inside/outside with and without tolerance
///  - inside/outside according to a given chi2 value
///
/// It also provides all the necessary tools for the individual implementations
/// in
/// the different SurfaceBounds classes.
///
/// @todo check if fast Sin/Cos ArcTan is necessary in field test
///
/// @todo Move the Covariance away from this into the method signature of the
/// call
/// @todo (short term) protect against lCovariance = nullptr acess

class BoundaryCheck
{
public:
  /// @brief nested enumerator for the boundary check Type
  enum BoundaryCheckType {
    absolute = 0,  ///< absolute check including tolerances
    chi2corr = 1   ///< relative (chi2 based) with full correlations
  };

  bool   checkLoc0;               ///< check local 1 coordinate
  bool   checkLoc1;               ///< check local 2 coordinate
  double toleranceLoc0;           ///< absolute tolerance in local 1 coordinate
  double toleranceLoc1;           ///< absolute tolerance in local 2 coordinate
  double nSigmas;                 ///< allowed sigmas for chi2 boundary check
  ActsSymMatrixD<2> lCovariance;  ///< local covariance matrix
  BoundaryCheckType bcType;       ///< the type how we check the boundary

  /// Constructor for single boolean behavious
  BoundaryCheck(bool sCheck);

  /// Constructor for tolerance based check
  /// @param chkL0 boolean directive to check first coordinate
  /// @param chkL1 boolean directive to check second coordinate
  /// @param tloc0 tolereance on the first parameter
  /// @param tloc1 tolereance on the second parameter
  BoundaryCheck(bool chkL0, bool chkL1, double tloc0 = 0., double tloc1 = 0.);

  /// Constructor for chi2 based check
  /// @param lCov the coverance matrix to be checked
  /// @param nsig number of sigma checked checked for the compatibility test
  /// @param chkL0 directive wheter the first coordinate is being checked
  /// @param chkL1 directive wheter the first coordinate is being checked
  BoundaryCheck(const ActsSymMatrixD<2>& lCov,
                double                   nsig  = 1.,
                bool                     chkL0 = true,
                bool                     chkL1 = true);

  /// Overwriting of the
  /// Conversion operator to bool
  operator bool() const { return (checkLoc0 || checkLoc1); }
  ///  Each Bounds has a method inside, which checks if a LocalPosition is
  ///  inside the bounds.
  ///  Inside can be called without/with boundary check */
  void
  ComputeKDOP(std::vector<Vector2D> v,
              std::vector<Vector2D> KDOPAxes,
              std::vector<KDOP>&    kdop) const;

  std::vector<Vector2D>
  EllipseToPoly(int resolution = 3) const;

  bool
  TestKDOPKDOP(std::vector<KDOP>& a, std::vector<KDOP>& b) const;
};

// does the conversion of an ellipse of height h and width w to an polygon with
// 4 + 4*resolution points
// @todo clean this code  & write documentation
inline std::vector<Vector2D>
BoundaryCheck::EllipseToPoly(int resolution) const
{
  const double h
      = (bcType == chi2corr) ? nSigmas * sqrt(lCovariance(1, 1)) : 0.;
  const double w
      = (bcType == chi2corr) ? nSigmas * sqrt(lCovariance(0, 0)) : 0.;

  // first add the four vertices
  std::vector<Vector2D> v((1 + resolution) * 4);
  Vector2D              p;
  p << w, 0;
  v.at(0) = p;
  p << -w, 0;
  v.at(1) = p;
  p << 0, h;
  v.at(2) = p;
  p << 0, -h;
  v.at(3) = p;

  // now add a number, equal to the resolution, of equidistant points  in each
  // quadrant
  // resolution == 3 seems to be a solid working point, but possibility open to
  // change to any number in the future
  Vector2D          t(1, 1);
  ActsSymMatrixD<2> t1;
  t1 << 1, 0, 0, -1;
  ActsSymMatrixD<2> t2;
  t2 << -1, 0, 0, -1;
  ActsSymMatrixD<2> t3;
  t3 << -1, 0, 0, 1;
  for (int i = 1; i <= resolution; i++) {
    double angle = M_PI_2 * i / (resolution + 1);
    t << w * std::sin(angle), h * std::cos(angle);
    v.at(i * 4 + 0) = t;
    v.at(i * 4 + 1) = t1 * t;
    v.at(i * 4 + 2) = t2 * t;
    v.at(i * 4 + 3) = t3 * t;
  }
  return v;
}

// calculates KDOP object from given polygon and set of axes
inline void
BoundaryCheck::ComputeKDOP(std::vector<Vector2D> v,
                           std::vector<Vector2D> KDOPAxes,
                           std::vector<KDOP>&    kdop) const
{
  // initialize KDOP to first point
  size_t k = KDOPAxes.size();
  for (size_t i = 0; i < k; i++) {
    kdop.at(i).max = KDOPAxes.at(i)(0, 0) * v.at(0)(0, 0)
        + KDOPAxes.at(i)(1, 0) * v.at(0)(1, 0);
    kdop.at(i).min = KDOPAxes.at(i)(0, 0) * v.at(0)(0, 0)
        + KDOPAxes.at(i)(1, 0) * v.at(0)(1, 0);
  }
  // now for each additional point, update KDOP bounds if necessary
  float value;
  for (size_t i = 1; i < v.size(); i++) {
    for (size_t j = 0; j < k; j++) {
      value = KDOPAxes.at(j)(0, 0) * v.at(i)(0, 0)
          + KDOPAxes.at(j)(1, 0) * v.at(i)(1, 0);
      if (value < kdop.at(j).min)
        kdop.at(j).min = value;
      else if (value > kdop.at(j).max)
        kdop.at(j).max = value;
    }
  }
}

// this is the method to actually check if two KDOPs overlap
inline bool
BoundaryCheck::TestKDOPKDOP(std::vector<KDOP>& a, std::vector<KDOP>& b) const
{
  size_t k = a.size();
  // check if any intervals are non-overlapping, return if so
  for (size_t i = 0; i < k; i++)
    if (a.at(i).min > b.at(i).max || a.at(i).max < b.at(i).min) return false;
  // all intervals are overlapping, so KDOPs must intersect
  return true;
}
}

#endif  // ACTS_SURFACES_BOUNDARYCHECK_H

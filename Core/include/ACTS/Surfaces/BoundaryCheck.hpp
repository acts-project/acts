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

#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/ParameterDefinitions.hpp"
// STD
#include <cmath>
#include <vector>

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

// struct needed for FastSinCos method (see below)
struct sincosCache
{
  double sinC;
  double cosC;
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
/// It also provides all the necessary tools for the individual implementations in
/// the different SurfaceBounds classes.
///
/// @TODO check if fast Sin/Cos ArcTan is necessary in field test
///
/// @TODO Move the Covariance away from this into the method signature of the call
 
class BoundaryCheck
{
  /// saves us a lot of function calls in the EllipseToPoly method
  /// @TODO fix it
    static constexpr double s_cos22 = 1.; // cos(22.5*M_PI/180.);
    static constexpr double s_cos45 = 1.; // cos(45.0*M_PI/180.);
    static constexpr double s_cos67 = 1.; // cos(67.5*M_PI/180.);

public:
    
  /// @enum nested enumerator for the boundary check Type
  enum BoundaryCheckType {
    absolute = 0,  ///< absolute check including tolerances
    chi2corr = 1   ///< relative (chi2 based) with full correlations
  };

  bool              checkLoc0;        ///< check local 1 coordinate
  bool              checkLoc1;        ///< check local 2 coordinate
  double            toleranceLoc0;    ///< absolute tolerance in local 1 coordinate
  double            toleranceLoc1;    ///< absolute tolerance in local 2 coordinate
  double            nSigmas;          ///< allowed sigmas for chi2 boundary check
  ActsSymMatrixD<2> lCovariance;      ///< local covariance matrix
  BoundaryCheckType bcType;           ///< the type how we check the boundary

  /// Constructor for single boolean behavious
  BoundaryCheck(bool sCheck)
    : checkLoc0(sCheck)
    , checkLoc1(sCheck)
    , toleranceLoc0(0.)
    , toleranceLoc1(0.)
    , nSigmas(-1)
    , lCovariance(ActsSymMatrixD<2>::Identity())
    , bcType(absolute)
  {}

  /// Constructor for tolerance based check
  /// @param chkL0 boolean directive to check first coordinate
  /// @param chkL1 boolean directive to check second coordinate
  /// @param tloc0 tolereance on the first parameter
  /// @param tloc1 tolereance on the second parameter
  BoundaryCheck(bool chkL0, bool chkL1, double tloc0 = 0., double tloc1 = 0.)
    : checkLoc0(chkL0)
    , checkLoc1(chkL1)
    , toleranceLoc0(tloc0)
    , toleranceLoc1(tloc1)
    , nSigmas(-1)
    , lCovariance(ActsSymMatrixD<2>::Identity())
    , bcType(absolute)
  {}

  /// Constructor for chi2 based check 
  /// @param lCov the coverance matrix to be checked 
  /// @param nsig number of sigma checked checked for the compatibility test
  /// @param chkL0 directive wheter the first coordinate is being checked
  /// @param chkL1 directive wheter the first coordinate is being checked
  BoundaryCheck(const ActsSymMatrixD<2>& lCov,
                double                   nsig  = 1.,
                bool                     chkL0 = true,
                bool                     chkL1 = true)
    : checkLoc0(chkL0)
    , checkLoc1(chkL1)
    , toleranceLoc0(0.)
    , toleranceLoc1(0.)
    , nSigmas(nsig)
    , lCovariance(lCov)
    , bcType(chi2corr)
  {}

  /// Overwriting of the 
  /// Conversion operator to bool
  operator bool() const { return (checkLoc0 || checkLoc1); }
  
  ///  Each Bounds has a method inside, which checks if a LocalPosition is inside the bounds.
  ///  Inside can be called without/with boundary check */
  void
  ComputeKDOP(std::vector<Vector2D> v,
              std::vector<Vector2D> KDOPAxes,
              std::vector<KDOP>&    kdop) const;

  std::vector<Vector2D>
  EllipseToPoly(int resolution = 3) const;

  bool
  TestKDOPKDOP(std::vector<KDOP>& a, std::vector<KDOP>& b) const;

  double
  FastArcTan(double x) const;

  sincosCache
  FastSinCos(double x) const;
};

/// should have maximum (average) error of 0.0015 (0.00089) radians or 0.0859
/// (0.0509) degrees, fine for us and much faster (>4 times)
inline double
BoundaryCheck::FastArcTan(double x) const
{
  double y;
  bool   complement = false;  // true if arg was >1
  bool   sign       = false;  // true if arg was < 0
  if (x < 0.) {
    x    = -x;
    sign = true;  // arctan(-x)=-arctan(x)
  }
  if (x > 1.) {
    x          = 1. / x;  // keep arg between 0 and 1
    complement = true;
  }
  y = M_PI_4 * x - x * (fabs(x) - 1) * (0.2447 + 0.0663 * fabs(x));
  if (complement) y = M_PI_2 - y;  // correct for 1/x if we did that
  if (sign) y       = -y;          // correct for negative arg
  return y;
}

/// should have maximum (average) error of 0.001 (0.0005) radians or 0.0573
/// (0.029) degrees, fine for us and much faster (>8 times)
inline sincosCache
BoundaryCheck::FastSinCos(double x) const
{
  sincosCache tmp;
  // always wrap input angle to -PI..PI
  if (x < -M_PI)
    x += 2. * M_PI;
  else if (x > M_PI)
    x -= 2. * M_PI;

  // compute sine
  if (x < 0.) {
    tmp.sinC = 1.27323954 * x + .405284735 * x * x;

    if (tmp.sinC < 0.)
      tmp.sinC = .225 * (tmp.sinC * -tmp.sinC - tmp.sinC) + tmp.sinC;
    else
      tmp.sinC = .225 * (tmp.sinC * tmp.sinC - tmp.sinC) + tmp.sinC;
  } else {
    tmp.sinC = 1.27323954 * x - 0.405284735 * x * x;

    if (tmp.sinC < 0.)
      tmp.sinC = .225 * (tmp.sinC * -tmp.sinC - tmp.sinC) + tmp.sinC;
    else
      tmp.sinC = .225 * (tmp.sinC * tmp.sinC - tmp.sinC) + tmp.sinC;
  }

  // compute cosine: sin(x + PI/2) = cos(x)
  x += M_PI_2;
  if (x > M_PI) x -= 2. * M_PI;

  if (x < 0.) {
    tmp.cosC = 1.27323954 * x + 0.405284735 * x * x;

    if (tmp.cosC < 0.)
      tmp.cosC = .225 * (tmp.cosC * -tmp.cosC - tmp.cosC) + tmp.cosC;
    else
      tmp.cosC = .225 * (tmp.cosC * tmp.cosC - tmp.cosC) + tmp.cosC;
  } else {
    tmp.cosC = 1.27323954 * x - 0.405284735 * x * x;

    if (tmp.cosC < 0.)
      tmp.cosC = .225 * (tmp.cosC * -tmp.cosC - tmp.cosC) + tmp.cosC;
    else
      tmp.cosC = .225 * (tmp.cosC * tmp.cosC - tmp.cosC) + tmp.cosC;
  }
  return tmp;
}

// does the conversion of an ellipse of height h and width w to an polygon with
// 4 + 4*resolution points
inline std::vector<Vector2D>
BoundaryCheck::EllipseToPoly(int resolution) const
{
  const double h = nSigmas * sqrt(lCovariance(1, 1));
  const double w = nSigmas * sqrt(lCovariance(0, 0));

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
  if (resolution != 3) {
    sincosCache scResult;
    for (int i = 1; i <= resolution; i++) {
      scResult = FastSinCos(M_PI_2 * i / (resolution + 1));
      t << w * scResult.sinC, h * scResult.cosC;
      v.at(i * 4 + 0) = t;
      v.at(i * 4 + 1) = t1 * t;
      v.at(i * 4 + 2) = t2 * t;
      v.at(i * 4 + 3) = t3 * t;
    }
  } else {
    t << w * s_cos22, h * s_cos67;
    v.at(4) = t;
    v.at(5) = t1 * t;
    v.at(6) = t2 * t;
    v.at(7) = t3 * t;
    t << w * s_cos45, h * s_cos45;
    v.at(8)  = t;
    v.at(9)  = t1 * t;
    v.at(10) = t2 * t;
    v.at(11) = t3 * t;
    t << w * s_cos67, h * s_cos22;
    v.at(12) = t;
    v.at(13) = t1 * t;
    v.at(14) = t2 * t;
    v.at(15) = t3 * t;
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
  unsigned int k = KDOPAxes.size();
  for (unsigned int i = 0; i < k; i++) {
    kdop.at(i).max = KDOPAxes.at(i)(0, 0) * v.at(0)(0, 0)
        + KDOPAxes.at(i)(1, 0) * v.at(0)(1, 0);
    kdop.at(i).min = KDOPAxes.at(i)(0, 0) * v.at(0)(0, 0)
        + KDOPAxes.at(i)(1, 0) * v.at(0)(1, 0);
  }
  // now for each additional point, update KDOP bounds if necessary
  float value;
  for (unsigned int i = 1; i < v.size(); i++) {
    for (unsigned int j = 0; j < k; j++) {
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
  int k = a.size();
  // check if any intervals are non-overlapping, return if so
  for (int i = 0; i < k; i++)
    if (a.at(i).min > b.at(i).max || a.at(i).max < b.at(i).min) return false;
  // all intervals are overlapping, so KDOPs must intersect
  return true;
}
}

#endif  // ACTS_SURFACES_BOUNDARYCHECK_H

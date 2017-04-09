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

  ActsSymMatrixD<2> lCovariance;  ///< local covariance matrix
  ActsSymMatrixD<2> m_weight;     ///< Weight matrix for metric
  double toleranceLoc0;           ///< absolute tolerance in local 1 coordinate
  double toleranceLoc1;           ///< absolute tolerance in local 2 coordinate
  double nSigmas;                 ///< allowed sigmas for chi2 boundary check
  bool   checkLoc0;               ///< check local 1 coordinate
  bool   checkLoc1;               ///< check local 2 coordinate
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
  /// @param lCov  the coverance matrix to be checked
  /// @param nsig  number of sigma checked checked for the compatibility test
  BoundaryCheck(const ActsSymMatrixD<2>& lCov, double nsig = 1);

  /// Return a new BoundaryCheck with updated covariance.
  /// @param jacobian Tranform Jacobian for the covariance
  /// @warning This currently only transforms the covariance and does not work
  ///          for the tolerance based check.
  BoundaryCheck
  transformed(const ActsMatrixD<2, 2>& jacobian) const;

  /// Check if the point is inside a polygon.
  ///
  /// @param point    Test point
  /// @param vertices Forward iterable container of convex polygon vertices
  ///
  /// The check takes into account whether tolerances or covariances are defined
  /// for the boundary check.
  template <typename Vector2DContainer>
  bool
  isInside(const Vector2D& point, const Vector2DContainer& vertices) const;
  /// Check if the point is inside a box aligned with the local axes.
  ///
  /// @param point   Test point
  /// @param loc0Min Lower bound along first axis
  /// @param loc0Max Upper bound along first axis
  /// @param loc1Min Lower bound along second axis
  /// @param loc1Max Upper bound along second axis
  ///
  /// The check takes into account whether tolerances or covariances are defined
  /// for the boundary check.
  bool
  isInside(const Vector2D& point,
           double          loc0Min,
           double          loc0Max,
           double          loc1Min,
           double          loc1Max) const;

  /// Calculate the signed, weighted, closest distance to a polygonal boundary.
  ///
  /// @param point Test point
  /// @param vertices Forward iterable container of convex polygon vertices
  /// @return Negative value if inside, positive if outside
  ///
  /// If a covariance is defined, the distance is the corresponding Mahalanobis
  /// distance. Otherwise, it is the Eucleadian distance.
  template <typename Vector2DContainer>
  double
  distance(const Vector2D& point, const Vector2DContainer& vertices) const;
  /// Calculate the signed, weighted, closest distance to an aligned box.
  ///
  /// @param point Test point
  /// @param loc0Min Minimal value along the first local axis
  /// @param loc0Max Maximal value along the first local axis
  /// @param loc1Min Minimal value along the first local axis
  /// @param loc1Max Maximal value along the first local axis
  /// @return Negative value if inside, positive if outside
  ///
  /// If a covariance is defined, the distance is the corresponding Mahalanobis
  /// distance. Otherwise, it is the Eucleadian distance.
  double
  distance(const Vector2D& point,
           double          loc0Min,
           double          loc0Max,
           double          loc1Min,
           double          loc1Max) const;

  /// Overwriting of the Conversion operator to bool
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

private:
  /// Check if the point is inside the polygon w/o any tolerances.
  template <typename Vector2DContainer>
  bool
  isInsidePolygon(const Vector2D&          point,
                  const Vector2DContainer& vertices) const;
  /// Check if the distance vector is within the absolute or relative limits.
  bool
  isTolerated(const Vector2D& delta) const;
  /// Compute vector norm based on the covariance.
  double
  squaredNorm(const Vector2D& x) const;
  /// Calculate the closest point on the polygon.
  template <typename Vector2DContainer>
  Vector2D
  computeClosestPointOnPolygon(const Vector2D&          point,
                               const Vector2DContainer& vertices) const;
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

inline Acts::BoundaryCheck::BoundaryCheck(bool sCheck)
  : lCovariance(ActsSymMatrixD<2>::Identity())
  , m_weight(ActsSymMatrixD<2>::Identity())
  , toleranceLoc0(0.)
  , toleranceLoc1(0.)
  , nSigmas(-1)
  , checkLoc0(sCheck)
  , checkLoc1(sCheck)
  , bcType(absolute)
{
}

inline Acts::BoundaryCheck::BoundaryCheck(bool   chkL0,
                                          bool   chkL1,
                                          double tloc0,
                                          double tloc1)
  : lCovariance(ActsSymMatrixD<2>::Identity())
  , m_weight(ActsSymMatrixD<2>::Identity())
  , toleranceLoc0(tloc0)
  , toleranceLoc1(tloc1)
  , nSigmas(-1)
  , checkLoc0(chkL0)
  , checkLoc1(chkL1)
  , bcType(absolute)
{
}

inline Acts::BoundaryCheck::BoundaryCheck(const ActsSymMatrixD<2>& lCov,
                                          double                   nsig)
  : lCovariance(lCov)
  , m_weight(lCov.inverse())
  , toleranceLoc0(0.)
  , toleranceLoc1(0.)
  , nSigmas(nsig)
  , checkLoc0(true)
  , checkLoc1(true)
  , bcType(chi2corr)
{
}

inline Acts::BoundaryCheck
Acts::BoundaryCheck::transformed(const ActsMatrixD<2, 2>& jacobian) const
{
  if (bcType == chi2corr) {
    return BoundaryCheck(jacobian * m_weight.inverse() * jacobian.transpose(),
                         nSigmas);
  } else {
    // project tolerances to the new system. depending on the jacobian we need
    // to check both tolerances, even when the initial check does not.
    Vector2D tol = jacobian * Vector2D(toleranceLoc0, toleranceLoc1);
    return BoundaryCheck(true, true, std::abs(tol[0]), std::abs(tol[1]));
  }
}

template <typename Vector2DContainer>
inline bool
Acts::BoundaryCheck::isInside(const Vector2D&          point,
                              const Vector2DContainer& vertices) const
{
  // a compatible point must be either completely in the polygon or on the
  // outside but within the defined tolerance relative to the closest point
  if (isInsidePolygon(point, vertices)) {
    return true;
  } else {
    auto closestPoint = computeClosestPointOnPolygon(point, vertices);
    return isTolerated(closestPoint - point);
  }
}

template <typename Vector2DContainer>
inline bool
Acts::BoundaryCheck::isInsidePolygon(const Vector2D&          point,
                                     const Vector2DContainer& vertices) const
{
  // when we move along the edges of a convex polygon, a point on the inside of
  // the polygon will always appear on the same side of each edge.
  // a point on the outside will switch sides at least once.

  // returns which side of the connecting line between `l0` and `l1` the point
  // `p` is on. computes the sign of the z-component of the cross-product
  // between the line normal vector and the vector from `l0` to `p`.
  auto lineSide = [&](auto&& l0, auto&& l1) {
    auto normal = l1 - l0;
    auto delta  = point - l0;
    return std::signbit((normal[0] * delta[1]) - (normal[1] * delta[0]));
  };

  auto     iv = std::begin(vertices);
  Vector2D l0 = *iv;
  Vector2D l1 = *(++iv);
  // use vertex0 to vertex1 to define reference sign and compare w/ all edges
  auto reference = lineSide(l0, l1);
  for (++iv; iv != std::end(vertices); ++iv) {
    l0 = l1;
    l1 = *iv;
    if (lineSide(l0, l1) != reference) {
      return false;
    }
  }
  // manual check for last edge from last vertex back to the first vertex
  if (lineSide(l1, *std::begin(vertices)) != reference) {
    return false;
  }
  // point was always on the same side. point must be inside.
  return true;
}

inline bool
Acts::BoundaryCheck::isInside(const Vector2D& point,
                              double          loc0Min,
                              double          loc0Max,
                              double          loc1Min,
                              double          loc1Max) const
{
  // TODO 2017-03-29 msmk: check open/close policy
  if ((loc0Min <= point[0]) && (point[0] < loc0Max) && (loc1Min <= point[1])
      && (point[1] < loc1Max)) {
    return true;
  } else {
    // TODO 2017-03-29 msmk: direct implementation for closest point on box
    Vector2D vertices[] = {{loc0Min, loc1Min},
                           {loc0Max, loc1Min},
                           {loc0Max, loc1Max},
                           {loc0Min, loc1Max}};
    Vector2D closestPoint = computeClosestPointOnPolygon(point, vertices);
    return isTolerated(closestPoint - point);
  }
}

inline bool
Acts::BoundaryCheck::isTolerated(const Vector2D& delta) const
{
  if (bcType == absolute) {
    bool insideLoc0 = !checkLoc0 || (std::abs(delta[0]) < toleranceLoc0);
    bool insideLoc1 = !checkLoc1 || (std::abs(delta[1]) < toleranceLoc1);
    return insideLoc0 && insideLoc1;
  } else {
    // 2d-Mahalanobis distance has an expectation value of 2
    return (squaredNorm(delta) < (2 * nSigmas));
  }
}

inline double
Acts::BoundaryCheck::squaredNorm(const Vector2D& x) const
{
  return (x.transpose() * lCovariance.inverse() * x).value();
}

template <typename Vector2DContainer>
inline double
Acts::BoundaryCheck::distance(const Acts::Vector2D&    point,
                              const Vector2DContainer& vertices) const
{
  // TODO 2017-04-06 msmk: this should be calculable directly
  double d = squaredNorm(point - computeClosestPointOnPolygon(point, vertices));
  d        = std::sqrt(d);
  return isInsidePolygonHard(point, vertices) ? -d : d;
}

inline double
Acts::BoundaryCheck::distance(const Acts::Vector2D& point,
                              double                loc0Min,
                              double                loc0Max,
                              double                loc1Min,
                              double                loc1Max) const
{
  // TODO 2017-04-04 msmk: direct implementation for closest point on box
  Vector2D vertices[] = {{loc0Min, loc1Min},
                         {loc0Max, loc1Min},
                         {loc0Max, loc1Max},
                         {loc0Min, loc1Max}};
  return distance(point, vertices);
}

template <typename Vector2DContainer>
inline Acts::Vector2D
Acts::BoundaryCheck::computeClosestPointOnPolygon(
    const Acts::Vector2D&    point,
    const Vector2DContainer& vertices) const
{
  // calculate the closest position on the segment between `l0` and `l1` to
  // the point as measured by the metric induced by the weight matrix
  auto closestOnSegment = [&](auto&& l0, auto&& l1) {
    // normal vector and position of the closest point along the normal
    auto n = l1 - l0;
    auto f = (n.transpose() * m_weight * n).value();
    auto u = -((l0 - point).transpose() * m_weight * n).value() / f;
    // u must be in [0, 1] to still be on the polygon segment
    return l0 + std::min(std::max(u, 0.0), 1.0) * n;
  };

  auto     iv      = std::begin(vertices);
  Vector2D l0      = *iv;
  Vector2D l1      = *(++iv);
  Vector2D closest = closestOnSegment(l0, l1);
  // Calculate the closest point on other connecting lines and compare distances
  for (++iv; iv != std::end(vertices); ++iv) {
    l0               = l1;
    l1               = *iv;
    Vector2D current = closestOnSegment(l0, l1);
    if (squaredNorm(current - point) < squaredNorm(closest - point)) {
      closest = current;
    }
  }
  // final edge from last vertex back to the first vertex
  Vector2D last = closestOnSegment(l1, *std::begin(vertices));
  if (squaredNorm(last - point) < squaredNorm(closest - point)) {
    closest = last;
  }
  return closest;
}

#endif  // ACTS_SURFACES_BOUNDARYCHECK_H

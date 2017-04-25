// This file is part of the ACTS project.
//
// Copyright (C) 2016-2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BoundaryCheck.hpp, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_BOUNDARYCHECK_H
#define ACTS_SURFACES_BOUNDARYCHECK_H

#include <cfloat>
#include <cmath>
#include <iterator>
#include <vector>

#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/ParameterDefinitions.hpp"

namespace Acts {

/// @class BoundaryCheck
///
/// The BoundaryCheck class provides boundary checks and distance calculations
/// for aligned box-like and polygonal boundaries on local surfaces.
/// Different types of boundary checks are supported and are transparently
/// selected when calling the `isInside(...)` and `distance(...)` methods:
///
/// -   Hard checks w/o any tolerances
/// -   Tolerance-based checks in one or in both local coordinates
/// -   Chi2-based checks based on a covariance matrix. Non-vanishing
///     correlations are correctly taken into account.
///
/// With a defined covariance matrix, the closest point and the distance are
/// not defined along the usual Euclidean metric, but by the Mahalanobis
/// distance induced by the the covariance.
class BoundaryCheck
{
public:
  /// Construct either hard cut in both dimensions or no cut at all.
  BoundaryCheck(bool check);
  /// Construct a tolerance based check.
  ///
  /// @param checkLocal0 Boolean directive to check coordinate 0
  /// @param checkLocal1 Boolean directive to check coordinate 1
  /// @param tolerance0 Tolerance along coordinate 0
  /// @param tolerance1 Tolerance along coordinate 1
  BoundaryCheck(bool   checkLocal0,
                bool   checkLocal1,
                double tolerance0 = 0,
                double tolerance1 = 0);
  /// Construct a chi2-based check.
  ///
  /// @param localCovariance Coverance matrix in local coordinates
  /// @param sigmaMax  Significance for the compatibility test
  BoundaryCheck(const ActsSymMatrixD<2>& localCovariance, double sigmaMax = 1);

  operator bool() const { return (m_type != Type::eNone); }
  bool operator!() const { return (m_type == Type::eNone); }

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

private:
  enum class Type {
    eNone,      ///< disable boundary check
    eAbsolute,  ///< absolute cut
    eChi2       ///< chi2-based cut with full correlations
  };

  /// Return a new BoundaryCheck with updated covariance.
  /// @param jacobian Tranform Jacobian for the covariance
  /// @warning This currently only transforms the covariance and does not work
  ///          for the tolerance based check.
  BoundaryCheck
  transformed(const ActsMatrixD<2, 2>& jacobian) const;

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

  /// metric weight matrix: identity for absolute mode or inverse covariance
  ActsSymMatrixD<2> m_weight;
  /// dual use: absolute tolerances or relative chi2/ sigma cut.
  Vector2D m_tolerance;
  Type     m_type;

  // To be able to use `transformed`
  friend class CylinderBounds;
  friend class DiscTrapezoidalBounds;
  // EllipseBounds needs a custom implementation
  friend class EllipseBounds;
};

}  // namespace Acts

inline Acts::BoundaryCheck::BoundaryCheck(bool check)
  : m_weight(ActsSymMatrixD<2>::Identity())
  , m_tolerance(0, 0)
  , m_type(check ? Type::eAbsolute : Type::eNone)
{
}

inline Acts::BoundaryCheck::BoundaryCheck(bool   checkLocal0,
                                          bool   checkLocal1,
                                          double tolerance0,
                                          double tolerance1)
  : m_weight(ActsSymMatrixD<2>::Identity())
  , m_tolerance(checkLocal0 ? tolerance0 : DBL_MAX,
                checkLocal1 ? tolerance1 : DBL_MAX)
  , m_type(Type::eAbsolute)
{
}

inline Acts::BoundaryCheck::BoundaryCheck(
    const ActsSymMatrixD<2>& localCovariance,
    double                   sigmaMax)
  : m_weight(localCovariance.inverse())
  , m_tolerance(sigmaMax, 0)
  , m_type(Type::eChi2)
{
}

inline Acts::BoundaryCheck
Acts::BoundaryCheck::transformed(const ActsMatrixD<2, 2>& jacobian) const
{
  BoundaryCheck bc = *this;
  if (m_type == Type::eAbsolute) {
    // project tolerances to the new system. depending on the jacobian we need
    // to check both tolerances, even when the initial check does not.
    bc.m_tolerance = (jacobian * m_tolerance).cwiseAbs();
  } else /* Type::eChi2 */ {
    bc.m_weight
        = (jacobian * m_weight.inverse() * jacobian.transpose()).inverse();
  }
  return bc;
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

template <typename Vector2DContainer>
inline double
Acts::BoundaryCheck::distance(const Acts::Vector2D&    point,
                              const Vector2DContainer& vertices) const
{
  // TODO 2017-04-06 msmk: this should be calculable directly
  double d = squaredNorm(point - computeClosestPointOnPolygon(point, vertices));
  d        = std::sqrt(d);
  return isInsidePolygon(point, vertices) ? -d : d;
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

inline bool
Acts::BoundaryCheck::isTolerated(const Vector2D& delta) const
{
  if (m_type == Type::eNone) {
    return true;
  } else if (m_type == Type::eAbsolute) {
    return (std::abs(delta[0]) <= m_tolerance[0])
        && (std::abs(delta[1]) <= m_tolerance[1]);
  } else /* Type::eChi2 */ {
    // Mahalanobis distances mean is 2 in 2-dim. cut is 1-d sigma.
    return (squaredNorm(delta) < (2 * m_tolerance[0]));
  }
}

inline double
Acts::BoundaryCheck::squaredNorm(const Vector2D& x) const
{
  return (x.transpose() * m_weight * x).value();
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

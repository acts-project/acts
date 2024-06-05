// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Surfaces/detail/VerticesHelper.hpp"

#include <cfloat>
#include <cmath>
#include <iterator>
#include <vector>

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
/// distance induced by the covariance.
class BoundaryCheck {
 public:
  /// Construct either hard cut in both dimensions or no cut at all.
  explicit BoundaryCheck(bool check);

  /// Construct a tolerance based check.
  ///
  /// @param checkLocal0 Boolean directive to check coordinate 0
  /// @param checkLocal1 Boolean directive to check coordinate 1
  /// @param tolerance0 Tolerance along coordinate 0
  /// @param tolerance1 Tolerance along coordinate 1
  BoundaryCheck(bool checkLocal0, bool checkLocal1, double tolerance0 = 0,
                double tolerance1 = 0);

  /// Construct a chi2-based check.
  ///
  /// @param localCovariance Coverance matrix in local coordinates
  /// @param sigmaMax  Significance for the compatibility test
  BoundaryCheck(const SquareMatrix2& localCovariance, double sigmaMax = 1);

  bool isEnabled() const { return m_type != Type::eNone; }

  /// Check if the point is inside a polygon.
  ///
  /// @param point    Test point
  /// @param vertices Forward iterable container of convex polygon vertices.
  ///                 Calling `std::begin`/ `std::end` on the container must
  ///                 return an iterator where `*it` must be convertible to
  ///                 an `Acts::Vector2`.
  ///
  /// The check takes into account whether tolerances or covariances are defined
  /// for the boundary check.
  template <typename Vector2Container>
  bool isInside(const Vector2& point, const Vector2Container& vertices) const;

  /// Check if the point is inside a box aligned with the local axes.
  ///
  /// @param point   Test point
  /// @param lowerLeft Minimal vertex of the box
  /// @param upperRight Maximal vertex of the box
  ///
  /// The check takes into account whether tolerances or covariances are defined
  /// for the boundary check.
  bool isInside(const Vector2& point, const Vector2& lowerLeft,
                const Vector2& upperRight) const;

  /// Calculate the signed, weighted, closest distance to a polygonal boundary.
  ///
  /// @param point Test point
  /// @param vertices Forward iterable container of convex polygon vertices.
  ///                 Calling `std::begin`/ `std::end` on the container must
  ///                 return an iterator where `*it` must be convertible to
  ///                 an `Acts::Vector2`.
  /// @return Negative value if inside, positive if outside
  ///
  /// If a covariance is defined, the distance is the corresponding Mahalanobis
  /// distance. Otherwise, it is the Eucleadian distance.
  template <typename Vector2Container>
  double distance(const Vector2& point, const Vector2Container& vertices) const;

  /// Calculate the signed, weighted, closest distance to an aligned box.
  ///
  /// @param point Test point
  /// @param lowerLeft Minimal vertex of the box
  /// @param upperRight Maximal vertex of the box
  /// @return Negative value if inside, positive if outside
  ///
  /// If a covariance is defined, the distance is the corresponding Mahalanobis
  /// distance. Otherwise, it is the Eucleadian distance.
  double distance(const Vector2& point, const Vector2& lowerLeft,
                  const Vector2& upperRight) const;

  enum class Type {
    eNone,      ///< disable boundary check
    eAbsolute,  ///< absolute cut
    eChi2       ///< chi2-based cut with full correlations
  };

  /// Broadcast the type
  Type type() const;

  // Broadcast the tolerance
  const Vector2& tolerance() const;

  // Return the covariance matrix
  SquareMatrix2 covariance() const;

 private:
  /// Return a new BoundaryCheck with updated covariance.
  /// @param jacobian Transform Jacobian for the covariance
  /// @warning This currently only transforms the covariance and does not work
  ///          for the tolerance based check.
  BoundaryCheck transformed(const ActsMatrix<2, 2>& jacobian) const;

  /// Check if the distance vector is within the absolute or relative limits.
  bool isTolerated(const Vector2& delta) const;

  /// Compute vector norm based on the covariance.
  double squaredNorm(const Vector2& x) const;

  /// Calculate the closest point on the polygon.
  template <typename Vector2Container>
  Vector2 computeClosestPointOnPolygon(const Vector2& point,
                                       const Vector2Container& vertices) const;

  /// metric weight matrix: identity for absolute mode or inverse covariance
  SquareMatrix2 m_weight;

  /// dual use: absolute tolerances or relative chi2/ sigma cut.
  Vector2 m_tolerance;
  Type m_type;

  // To access the m_type
  friend class CylinderBounds;
  friend class RectangleBounds;
  // To be able to use `transformed`
  friend class CylinderBounds;
  friend class DiscTrapezoidBounds;
  // EllipseBounds needs a custom implementation
  friend class EllipseBounds;
};

}  // namespace Acts

inline Acts::BoundaryCheck::Type Acts::BoundaryCheck::type() const {
  return m_type;
}

inline const Acts::Vector2& Acts::BoundaryCheck::tolerance() const {
  return m_tolerance;
}

inline Acts::SquareMatrix2 Acts::BoundaryCheck::covariance() const {
  return m_weight.inverse();
}

inline Acts::BoundaryCheck::BoundaryCheck(bool check)
    : m_weight(SquareMatrix2::Identity()),
      m_tolerance(0, 0),
      m_type(check ? Type::eAbsolute : Type::eNone) {}

inline Acts::BoundaryCheck::BoundaryCheck(bool checkLocal0, bool checkLocal1,
                                          double tolerance0, double tolerance1)
    : m_weight(SquareMatrix2::Identity()),
      m_tolerance(checkLocal0 ? tolerance0 : DBL_MAX,
                  checkLocal1 ? tolerance1 : DBL_MAX),
      m_type(Type::eAbsolute) {}

inline Acts::BoundaryCheck::BoundaryCheck(const SquareMatrix2& localCovariance,
                                          double sigmaMax)
    : m_weight(localCovariance.inverse()),
      m_tolerance(sigmaMax, 0),
      m_type(Type::eChi2) {}

inline Acts::BoundaryCheck Acts::BoundaryCheck::transformed(
    const ActsMatrix<2, 2>& jacobian) const {
  BoundaryCheck bc = *this;
  if (m_type == Type::eAbsolute) {
    // project tolerances to the new system. depending on the jacobian we need
    // to check both tolerances, even when the initial check does not.
    bc.m_tolerance = (jacobian * m_tolerance).cwiseAbs();
  } else /* Type::eChi2 */ {
    bc.m_weight =
        (jacobian * m_weight.inverse() * jacobian.transpose()).inverse();
  }
  return bc;
}

template <typename Vector2Container>
inline bool Acts::BoundaryCheck::isInside(
    const Vector2& point, const Vector2Container& vertices) const {
  if (m_type == Type::eNone) {
    // The null boundary check always succeeds
    return true;
  } else if (detail::VerticesHelper::isInsidePolygon(point, vertices)) {
    // If the point falls inside the polygon, the check always succeeds
    return true;
  } else if (m_tolerance == Vector2(0., 0.)) {
    // Outside of the polygon, since we've eliminated the case of an absence of
    // check above, we know we'll always fail if the tolerance is zero.
    //
    // This allows us to avoid the expensive computeClosestPointOnPolygon
    // computation in this simple case.
    //
    // TODO: When tolerance is not 0, we could also avoid this computation in
    //       some cases by testing against a bounding box of the polygon, padded
    //       on each side with our tolerance. Check if this optimization is
    //       worthwhile in some production workflows, and if so implement it.
    return false;
  } else {
    // We are outside of the polygon, but there is a tolerance. Must find what
    // the closest point on the polygon is and check if it's within tolerance.
    auto closestPoint = computeClosestPointOnPolygon(point, vertices);
    return isTolerated(closestPoint - point);
  }
}

inline bool Acts::BoundaryCheck::isInside(const Vector2& point,
                                          const Vector2& lowerLeft,
                                          const Vector2& upperRight) const {
  if (detail::VerticesHelper::isInsideRectangle(point, lowerLeft, upperRight)) {
    return true;
  } else {
    Vector2 closestPoint;

    if (m_type == Type::eNone || m_type == Type::eAbsolute) {
      // absolute, can calculate directly
      closestPoint =
          detail::VerticesHelper::computeEuclideanClosestPointOnRectangle(
              point, lowerLeft, upperRight);

    } else /* Type::eChi2 */ {
      // need to calculate by projection and squarednorm
      Vector2 vertices[] = {{lowerLeft[0], lowerLeft[1]},
                            {upperRight[0], lowerLeft[1]},
                            {upperRight[0], upperRight[1]},
                            {lowerLeft[0], upperRight[1]}};
      closestPoint = computeClosestPointOnPolygon(point, vertices);
    }

    return isTolerated(closestPoint - point);
  }
}

template <typename Vector2Container>
inline double Acts::BoundaryCheck::distance(
    const Acts::Vector2& point, const Vector2Container& vertices) const {
  // TODO 2017-04-06 msmk: this should be calculable directly
  double d = squaredNorm(point - computeClosestPointOnPolygon(point, vertices));
  d = std::sqrt(d);
  return detail::VerticesHelper::isInsidePolygon(point, vertices) ? -d : d;
}

inline double Acts::BoundaryCheck::distance(const Acts::Vector2& point,
                                            const Vector2& lowerLeft,
                                            const Vector2& upperRight) const {
  if (m_type == Type::eNone || m_type == Type::eAbsolute) {
    // compute closest point on box
    double d = (point -
                detail::VerticesHelper::computeEuclideanClosestPointOnRectangle(
                    point, lowerLeft, upperRight))
                   .norm();
    return detail::VerticesHelper::isInsideRectangle(point, lowerLeft,
                                                     upperRight)
               ? -d
               : d;

  } else /* Type::eChi2 */ {
    Vector2 vertices[] = {{lowerLeft[0], lowerLeft[1]},
                          {upperRight[0], lowerLeft[1]},
                          {upperRight[0], upperRight[1]},
                          {lowerLeft[0], upperRight[1]}};
    return distance(point, vertices);
  }
}

inline bool Acts::BoundaryCheck::isTolerated(const Vector2& delta) const {
  if (m_type == Type::eNone) {
    return true;
  } else if (m_type == Type::eAbsolute) {
    return (std::abs(delta[0]) <= m_tolerance[0]) &&
           (std::abs(delta[1]) <= m_tolerance[1]);
  } else /* Type::eChi2 */ {
    // Mahalanobis distances mean is 2 in 2-dim. cut is 1-d sigma.
    return (squaredNorm(delta) < (2 * m_tolerance[0]));
  }
}

inline double Acts::BoundaryCheck::squaredNorm(const Vector2& x) const {
  return (x.transpose() * m_weight * x).value();
}

template <typename Vector2Container>
inline Acts::Vector2 Acts::BoundaryCheck::computeClosestPointOnPolygon(
    const Acts::Vector2& point, const Vector2Container& vertices) const {
  return detail::VerticesHelper::computeClosestPointOnPolygon(point, vertices,
                                                              m_weight);
}

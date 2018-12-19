// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// Intersection.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <limits>
#include "Definitions.hpp"

namespace Acts {

/// A function typedef for the intersection correction
using CorrFnc = std::function<bool(Vector3D&, Vector3D&, double&)>;

///  @struct Intersection
///
///  intersection struct used for position
struct Intersection
{
  Vector3D position;      ///< position of the intersection
  double   pathLength;    ///< path length to the intersection (if valid)
  double   distance{0.};  ///< remaining distance (if not valid)
  bool     valid{false};  ///< validiaty boolean

  /// Constructor with arguments
  ///
  /// @param sinter is the position of the intersection
  /// @param slength is the path length to the intersection
  /// @param svalid is a boolean indicating if intersection is valid
  /// @param dist is the distance to the closes surface boundary
  Intersection(const Vector3D& sinter,
               double          slength,
               bool            svalid,
               double          dist = 0.)
    : position(sinter), pathLength(slength), distance(dist), valid(svalid)
  {
  }

  /// Default constructor
  Intersection()
    : position(Vector3D(0., 0., 0.))
    , pathLength(std::numeric_limits<double>::infinity())
  {
  }

  /// Bool() operator for validity checking
  explicit operator bool() const { return valid; }

  /// Smaller operator for sorting,
  /// - it respects the validity of the intersection
  /// @param si is the intersection for testing
  bool
  operator<(const Intersection& si) const
  {
    if (!valid) {
      return false;
    }
    // now check the pathLength
    if (si.valid) {
      return (pathLength < si.pathLength);
    }
    // the current path length wins
    return true;
  }

  /// Greater operator for sorting,
  /// - it respects the validity of the intersection
  /// @param si is the intersection for testing
  bool
  operator>(const Intersection& si) const
  {
    if (!valid) {
      return false;
    }
    // now check the pathLength
    if (si.valid) {
      return (pathLength > si.pathLength);
    }
    // the current path length wins
    return true;
  }
};

/// class extensions to return also the object - can be merged with
/// FullIntersection
template <typename object_t>
class ObjectIntersection
{
public:
  Intersection    intersection{};   ///< the intersection iself
  const object_t* object{nullptr};  ///< the object that was intersected
  const object_t* representation{
      nullptr};  ///< the representation of the object
  NavigationDirection pDirection{
      anyDirection};  ///< the direction in which it was taken

  /// Default constructor
  ObjectIntersection() = default;

  /// Object intersection
  ///
  /// @param sInter is the intersection
  /// @param sObject is the object to be instersected
  /// @param dir is the direction of the intersection
  ObjectIntersection(const Intersection& sInter,
                     const object_t*     sObject,
                     NavigationDirection dir = forward)
    : intersection(sInter)
    , object(sObject)
    , representation(sObject)
    , pDirection(dir)
  {
  }

  /// Bool() operator for validity checking
  explicit operator bool() const { return intersection.valid; }

  /// @brief smaller operator for ordering & sorting
  ///
  /// @param oi is the source intersection for comparison
  bool
  operator<(const ObjectIntersection<object_t>& oi) const
  {
    return (intersection < oi.intersection);
  }

  /// @brief greater operator for ordering & sorting
  ///
  /// @param oi is the source intersection for comparison
  bool
  operator>(const ObjectIntersection<object_t>& oi) const
  {
    return (intersection > oi.intersection);
  }
};

/// @brief Class with full intersection information
///
/// It contains the interscetion, the object that was intersected
/// and the representation of the object (identical if doesn't differ)
///
/// @tparam object_t Type of the object to be intersected
/// @tparam representation_t Type of the representation
template <typename object_t, typename representation_t>
class FullIntersection
{
public:
  Intersection            intersection;    ///< the intersection iself
  const object_t*         object;          ///< the object that was intersected
  const representation_t* representation;  ///< the represenation of the object
  NavigationDirection     pDirection;  ///< the direction in which it was taken

  /// Full intersection constructor
  ///
  /// @param sInter is the intersection struct
  /// @param sObject is the intersected object
  /// @param sRepresentation is the surface representation of the object
  /// @param dir is the direction
  ///
  FullIntersection(const Intersection&     sInter,
                   const object_t*         sObject,
                   const representation_t* sRepresentation,
                   NavigationDirection     dir = forward)
    : intersection(sInter)
    , object(sObject)
    , representation(sRepresentation)
    , pDirection(dir)
  {
  }

  /// Bool() operator for validity checking
  explicit operator bool() const { return intersection.valid; }

  /// @brief smaller operator for ordering & sorting
  ///
  /// @param fi is the full intersection to be tested
  bool
  operator<(const FullIntersection<object_t, representation_t>& fi) const
  {
    return (intersection < fi.intersection);
  }

  /// @brief greater operator for ordering & sorting
  ///
  /// @param fi is the full intersection to be tested
  bool
  operator>(const FullIntersection<object_t, representation_t>& fi) const
  {
    return (intersection > fi.intersection);
  }
};

struct SameSurfaceIntersection
{
  /// @brief comparison operator
  ///
  /// This is a struct to pick out intersection with identical surfaces
  ///
  /// @tparam intersection_t Type of the intersection object
  /// @param i1 First intersection to test
  /// @param i2 Second intersection to test
  template <typename intersection_t>
  bool
  operator()(const intersection_t& i1, const intersection_t& i2) const
  {
    return (i1.object == i2.object);
  }
};

/// @brief Void Direction corrector
///
/// This is used to evaluate a modified
/// intersection (e.g. curvature updated)
struct VoidIntersectionCorrector
{

  // Void Corrector default constructor
  VoidIntersectionCorrector() = default;

  // Void Corrector parameter constructor
  VoidIntersectionCorrector(const Vector3D& /*unused*/,
                            const Vector3D& /*unused*/,
                            double /*unused*/)
  {
  }

  /// Boolean() operator - returns false for void modifier
  explicit operator bool() const { return false; }

  /// empty correction interface
  bool
  operator()(Vector3D& /*unused*/, Vector3D& /*unused*/, double /*unused*/)
  {
    return false;
  }

  /// Step modification call
  ///
  /// @stay put and don't do antyhing
  template <typename step_t>
  bool
  operator()(step_t& /*unused*/) const
  {
    return false;
  }
};
}

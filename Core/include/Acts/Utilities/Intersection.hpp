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
  Vector3D position;    ///< position of the intersection
  double   pathLength;  ///< path length to the intersection (if valid)
  double   distance;    ///< remaining distance (if not valid)
  bool     valid;       ///< validiaty boolean

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
    , distance(0.)
    , valid(false)
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
    if (!valid) return false;
    // now check the pathLenght
    if (si.valid) return (pathLength < si.pathLength);
    // the current path length wins
    return true;
  }

  /// Greater operator for sorting,
  /// - it respects the validity of the intersection
  /// @param si is the intersection for testing
  bool
  operator>(const Intersection& si) const
  {
    if (!valid) return false;
    // now check the pathLenght
    if (si.valid) return (pathLength > si.pathLength);
    // the current path length wins
    return true;
  }
};

/// class extensions to return also the object
template <typename object_t>
class ObjectIntersection
{
public:
  Intersection        intersection;  ///< the intersection iself
  const object_t*     object;        ///< the object that was intersected
  NavigationDirection pDirection;    ///< the direction in which it was taken

  /// Default constructor
  ObjectIntersection()
    : intersection(), object(nullptr), pDirection(anyDirection)
  {
  }

  /// Object intersection
  ///
  /// @param sInter is the intersection
  /// @param sObject is the object to be instersected
  /// @param dir is the direction of the intersection
  ObjectIntersection(const Intersection& sInter,
                     const object_t*     sObject,
                     NavigationDirection dir = forward)
    : intersection(sInter), object(sObject), pDirection(dir)
  {
  }

  /// Bool() operator for validity checking
  explicit operator bool() const { return intersection.valid; }

  /// smaller operator for ordering & sorting
  ///
  /// @param oi is the source intersection for comparison
  bool
  operator<(const ObjectIntersection<object_t>& oi) const
  {
    return (intersection < oi.intersection);
  }

  /// greater operator for ordering & sorting
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

  /// Smaller operator for ordering & sorting
  ///
  /// @param fi is the full intersection to be tested
  bool
  operator<(const FullIntersection<object_t, representation_t>& fi) const
  {
    return (intersection < fi.intersection);
  }

  /// Greater operator for ordering & sorting
  ///
  /// @param fi is the full intersection to be tested
  bool
  operator>(const FullIntersection<object_t, representation_t>& fi) const
  {
    return (intersection > fi.intersection);
  }
};

/// @brief Void Direction corrector
///
/// This is used to evaluate a modified
/// intersection (e.g. curvature updated)
struct VoidCorrector
{

  // Void Corrector default constructor
  VoidCorrector() {}

  // Void Corrector parameter constructor
  VoidCorrector(const Vector3D&, const Vector3D&, double) {}

  /// Boolean() operator - returns false for void modifier
  explicit operator bool() const { return false; }

  /// empty correction interface
  bool
  operator()(Vector3D&, Vector3D&, double)
  {
    return false;
  }
};
}

// This file is part of the Acts project.
//
// Copyright (C) 2016-2017 Acts project team
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

///  @struct Intersection
///
///  intersection struct used for position
struct Intersection
{
  Vector3D position;  // position of the intersection
  double   pathLength;
  double   distance;
  bool     valid;

  /// Constructor with argoments
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

  /// Smaller operator for sorting
  ///
  /// @param si is the intersection for testing
  bool
  operator<(const Intersection& si) const
  {
    return (valid && pathLength < si.pathLength);
  }
};

/// class extensions to return also the object
template <class T>
class ObjectIntersection
{
public:
  Intersection intersection;
  const T*     object;
  int          pDirection;

  /// Default constructor
  ObjectIntersection() : intersection(), object(nullptr), pDirection(0) {}
  /// Object intersection
  ///
  /// @param sInter is the intersection
  /// @param sObject is the object to be instersected
  /// @param dir is the direction of the intersection
  ObjectIntersection(const Intersection& sInter, const T* sObject, int dir = 1)
    : intersection(sInter), object(sObject), pDirection(dir)
  {
  }

  /// Bool() operator for validity checking
  explicit operator bool() const { return intersection.valid; }

  /// smaller operator for ordering & sorting
  ///
  /// @param oi is the source intersection for comparison
  bool
  operator<(const ObjectIntersection<T>& oi) const
  {
    return (intersection < oi.intersection);
  }
};

/// Class extension to return the object, a represenation & the result
template <class T, class R, class S>
class FullIntersection
{
public:
  Intersection intersection;
  const T*     object;
  const R*     representation;
  const S*     result;  //@todo: remove bare pointer
  int          pDirection;

  /// Full intersection constructor
  ///
  /// @param sInter is the intersection struct
  /// @param sObject is the intersected object
  /// @param sRepresentation is the surface representation of the object
  /// @param sResult is the type of result: neutral, charged TP e.g.
  /// @param dir is the direction
  ///
  /// @todo use unique_ptr for result !
  FullIntersection(const Intersection& sInter,
                   const T*            sObject,
                   const R*            sRepresentation,
                   const S*            sResult = nullptr,
                   int                 dir     = 1)
    : intersection(sInter)
    , object(sObject)
    , representation(sRepresentation)
    , result(sResult)
    , pDirection(dir)
  {
  }

  /// Bool() operator for validity checking
  explicit operator bool() const { return intersection.valid; }

  /// Smaller operator for ordering & sorting
  ///
  /// @param fi is the full intersection to be tested
  bool
  operator<(const FullIntersection<T, R, S>& fi) const
  {
    return (intersection < fi.intersection);
  }
};
}
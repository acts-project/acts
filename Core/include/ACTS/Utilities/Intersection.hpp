// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// Intersection.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYUTILS_INTERSECTION_H
#define ACTS_GEOMETRYUTILS_INTERSECTION_H 1

#include "Definitions.hpp"

namespace Acts {

///  @struct Intersection
///
///  intersection struct used for position
struct Intersection
{
  Vector3D position;
  double   pathLength;
  double   distance;
  bool     valid;

  /// Constructor with argoments
  ///
  /// @param sinter is the position of the intersection
  /// @param slength is the path length to the intersection
  /// @param svalid is a boolean indicating if intersection is valid
  /// @param dist is the distance to the surface, e.g. when outside bounds
  Intersection(const Vector3D& sinter,
               double          slength,
               bool            svalid,
               double          dist = 0.)
    : position(sinter), pathLength(slength), distance(dist), valid(svalid)
  {
  }

  Intersection()
    : position(Vector3D(0., 0., 0.)), pathLength(0.), distance(0.), valid(false)
  {
  }

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
  Intersection     intersection;
  mutable const T* object;
  int              pDirection;

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
  Intersection     intersection;
  mutable const T* object;
  mutable const R* representation;
  mutable const S* result;
  int              pDirection;

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
                   const S*            sResult,
                   int                 dir = 1)
    : intersection(sInter)
    , object(sObject)
    , representation(sRepresentation)
    , result(sResult)
    , pDirection(dir)
  {
  }

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

#endif  // ACTS_GEOMETRYUTILS_INTERSECTION_H

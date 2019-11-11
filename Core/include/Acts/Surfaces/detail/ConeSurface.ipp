// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ConeSurface.ipp, Acts project
///////////////////////////////////////////////////////////////////

inline Intersection ConeSurface::intersectionEstimate(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& direction, const BoundaryCheck& bcheck) const {
  // check if you need
  bool needsTransform = (Surface::m_transform) ? true : false;

  // Create the points
  Vector3D point1 = position;
  Vector3D dir1 = direction;

  // What you need at the and
  Vector3D solution(0, 0, 0);
  double path = 0.;
  Intersection::Status status = Intersection::Status::unreachable;

  if (needsTransform) {
    Transform3D invTrans = transform(gctx).inverse();
    point1 = invTrans * position;
    dir1 = invTrans.linear() * direction;
  }

  // See file header for the formula derivation
  double tan2Alpha = bounds().tanAlpha() * bounds().tanAlpha(),
         A = dir1.x() * dir1.x() + dir1.y() * dir1.y() -
             tan2Alpha * dir1.z() * dir1.z(),
         B = 2 * (dir1.x() * point1.x() + dir1.y() * point1.y() -
                  tan2Alpha * dir1.z() * point1.z()),
         C = point1.x() * point1.x() + point1.y() * point1.y() -
             tan2Alpha * point1.z() * point1.z();
  if (A == 0.) {
    A += 1e-16;  // avoid division by zero
  }

  detail::RealQuadraticEquation solns(A, B, C);

  if (solns.solutions != 0) {
    double t1 = solns.first;
    Vector3D soln1Loc(point1 + t1 * dir1);
    // set the validity flag
    status = Intersection::Status::reachable;
    // there's only one solution for this
    if (solns.solutions == 1) {
      // set the solution
      solution = soln1Loc;
      path = t1;
      // check the validity given the navigation direction
      status = Intersection::Status::reachable;
      // there's two solutions
    } else if (solns.solutions == 2) {
      // get the second solution
      double t2 = solns.second;
      Vector3D soln2Loc(point1 + t2 * direction);
      // both solutions have the same sign - or you don't care
      // @TODO needs two solution
      if (t1 * t1 < t2 * t2) {
        solution = soln1Loc;
        path = t1;
      } else {
        solution = soln2Loc;
        path = t2;
      }
    }
  }

  // Transform back if needed
  if (m_transform) {
    solution = transform(gctx) * solution;
  }
  // Check validity
  if (status != Intersection::Status::unreachable and bcheck and
      not isOnSurface(gctx, solution, direction, bcheck)) {
    status = Intersection::Status::missed;
  }

  // Set the result navigation direction
  return Intersection(solution, path, status);
}

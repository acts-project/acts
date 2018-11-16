// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ConeSurface.ipp, Acts project
///////////////////////////////////////////////////////////////////

inline Intersection
ConeSurface::intersectionEstimate(const Vector3D&      gpos,
                                  const Vector3D&      gmom,
                                  NavigationDirection  navDir,
                                  const BoundaryCheck& bcheck,
                                  CorrFnc              correct) const
{

  // check if you need
  bool needsTransform = (Surface::m_transform) ? true : false;

  // create the points
  Vector3D point1     = gpos;
  Vector3D initialdir = gmom.normalized();
  Vector3D direction  = initialdir;

  // what you need at the and
  Vector3D solution(0, 0, 0);
  double   path  = 0.;
  bool     valid = false;

  // break condition for the loop
  bool correctionDone = false;

  // make the loop including the potential correction
  do {

    if (needsTransform) {
      Transform3D invTrans = transform().inverse();
      point1               = invTrans * gpos;
      direction            = invTrans.linear() * direction;
    }

    // see the header for the formula derivation
    double tan2Alpha = bounds().tanAlpha() * bounds().tanAlpha(),
           A = direction.x() * direction.x() + direction.y() * direction.y()
        - tan2Alpha * direction.z() * direction.z(),
           B = 2 * (direction.x() * point1.x() + direction.y() * point1.y()
                    - tan2Alpha * direction.z() * point1.z()),
           C = point1.x() * point1.x() + point1.y() * point1.y()
        - tan2Alpha * point1.z() * point1.z();
    if (A == 0.) {
      A += 1e-16;  // avoid division by zero
    }

    // use Andreas' quad solver, much more stable than what I wrote
    detail::RealQuadraticEquation solns(A, B, C);

    if (solns.solutions != 0) {
      double   t1 = solns.first;
      Vector3D soln1Loc(point1 + t1 * direction);

      // there's only one solution for this
      if (solns.solutions == 1) {
        // set the solution
        solution = soln1Loc;
        path     = t1;
        // check the validity given the navigation direction
        valid = (t1 * navDir >= 0.);
        // there's two solutions
      } else if (solns.solutions == 2) {
        // get the second solution
        double   t2 = solns.second;
        Vector3D soln2Loc(point1 + t2 * direction);
        // both solutions have the same sign - or you don't care
        // then take the closer one
        if (t1 * t2 > 0. || (navDir == 0)) {
          if (t1 * t1 < t2 * t2) {
            solution = soln1Loc;
            path     = t1;
          } else {
            solution = soln2Loc;
            path     = t2;
          }
        } else {
          if (navDir * t1 > 0.) {
            solution = soln1Loc;
            path     = t1;
          } else {
            solution = soln2Loc;
            path     = t2;
          }
        }
      }
    }
    // set the validity flag
    valid = (navDir * path >= 0.);

    // break if there's nothing to correct
    if (correct && !correctionDone) {
      // reset to initial position and direction
      point1    = gpos;
      direction = initialdir;
      if (correct(point1, direction, path)) {
        correctionDone = true;
      } else {
        break;
      }
    } else {
      break;
    }
  } while (true);

  // transform back if needed
  if (m_transform) {
    solution = transform() * solution;
  }
  // check validity
  valid = bcheck ? (valid && isOnSurface(solution, gmom, bcheck)) : valid;
  // set the result navigation direction
  return Intersection(solution, path, valid);
}

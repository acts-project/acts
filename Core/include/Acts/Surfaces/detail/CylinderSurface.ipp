// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderSurface.ipp, Acts project
///////////////////////////////////////////////////////////////////

inline const Vector3D
CylinderSurface::rotSymmetryAxis() const
{
  // fast access via tranform matrix (and not rotation())
  return transform().matrix().block<3, 1>(0, 2);
}

inline Intersection
CylinderSurface::intersectionEstimate(const Vector3D&      gpos,
                                      const Vector3D&      gdir,
                                      NavigationDirection  navDir,
                                      const BoundaryCheck& bcheck,
                                      CorrFnc              correct) const
{

  // create line parameters
  Vector3D lpos = gpos;
  Vector3D ldir = gdir;
  // minimize the call to transform()
  const auto& tMatrix = transform().matrix();
  Vector3D    caxis   = tMatrix.block<3, 1>(0, 2).transpose();
  Vector3D    ccenter = tMatrix.block<3, 1>(0, 3).transpose();
  // what you need at the and
  Vector3D solution(0, 0, 0);
  double   path = 0.;

  // lemma : the solver -> should catch current values
  auto solve = [&solution, &path, &lpos, &ldir, &ccenter, &caxis, &navDir](
      double R) -> bool {
    // check documentation for explanation
    Vector3D pc    = lpos - ccenter;
    Vector3D pcXcd = pc.cross(caxis);
    Vector3D ldXcd = ldir.cross(caxis);
    double   a     = ldXcd.dot(ldXcd);
    double   b     = 2. * (ldXcd.dot(pcXcd));
    double   c     = pcXcd.dot(pcXcd) - (R * R);
    // and solve the qaudratic equation - todo, validity check
    detail::RealQuadraticEquation qe(a, b, c);
    // check how many solution you have
    if (qe.solutions == 0) {
      return false;
    }
    // chose the solution
    path = ((navDir == 0) || qe.first * qe.second > 0.)
        ? (qe.first * qe.first < qe.second * qe.second ? qe.first : qe.second)
        : (navDir * qe.first >= 0. ? qe.first : qe.second);
    // return the solution
    solution = lpos + path * ldir;
    // is valid if it goes into the right direction
    return (path * navDir >= 0.);
  };

  // solve for radius R
  double R     = bounds().r();
  bool   valid = solve(R);
  // if configured, correct and solve again
  if (correct && correct(lpos, ldir, path)) {
    valid = solve(R);
  }
  // update for inside if requested :
  // @todo fix this : fast inside bounds check needed
  valid = bcheck ? (valid && isOnSurface(solution, gdir, bcheck)) : valid;
  // now return
  return Intersection(solution, path, valid);
}

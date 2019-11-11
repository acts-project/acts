// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderSurface.ipp, Acts project
///////////////////////////////////////////////////////////////////

inline const Vector3D CylinderSurface::rotSymmetryAxis(
    const GeometryContext& gctx) const {
  // fast access via tranform matrix (and not rotation())
  return transform(gctx).matrix().block<3, 1>(0, 2);
}

inline Intersection CylinderSurface::intersectionEstimate(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& direction, const BoundaryCheck& bcheck) const {
  // create line parameters
  Vector3D lpos = position;
  Vector3D ldir = direction;
  // minimize the call to transform()
  const auto& tMatrix = transform(gctx).matrix();
  Vector3D caxis = tMatrix.block<3, 1>(0, 2).transpose();
  Vector3D ccenter = tMatrix.block<3, 1>(0, 3).transpose();
  // what you need at the and
  Vector3D solution(0, 0, 0);
  double path = 0.;

  // lemma : the solver ----- encapsulated
  auto solve = [&solution, &path, &lpos, &ldir, &ccenter,
                &caxis](double R) -> Intersection::Status {
    // Check documentation for explanation
    Vector3D pc = lpos - ccenter;
    Vector3D pcXcd = pc.cross(caxis);
    Vector3D ldXcd = ldir.cross(caxis);
    double a = ldXcd.dot(ldXcd);
    double b = 2. * (ldXcd.dot(pcXcd));
    double c = pcXcd.dot(pcXcd) - (R * R);
    // and solve the qaudratic equation - todo, validity check
    detail::RealQuadraticEquation qe(a, b, c);
    // check how many solution you have
    if (qe.solutions == 0) {
      return Intersection::Status::unreachable;
    }
    // @TODO this needs some thinking ... we need to provide both solutions
    // Chose the solution
    path = qe.first * qe.first < qe.second * qe.second ? qe.first : qe.second;
    // return the solution
    solution = lpos + path * ldir;
    // is valid if it goes into the right direction
    return Intersection::Status::reachable;
  };
  // ------

  // Solve for radius R
  double R = bounds().r();
  Intersection::Status status = solve(R);
  // Boundary check necessary
  if (status != Intersection::Status::unreachable and bcheck and
      not isOnSurface(gctx, solution, direction, bcheck)) {
    status = Intersection::Status::missed;
  }
  // Now return the solution
  return Intersection(solution, path, status);
}

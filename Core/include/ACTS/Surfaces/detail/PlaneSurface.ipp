// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// PlaneSurface.ipp, ACTS project
///////////////////////////////////////////////////////////////////

inline const Vector3D
PlaneSurface::normal(const Vector2D&) const
{
  // fast access via tranform matrix (and not rotation())
  const auto& tMatrix = transform().matrix();
  return Vector3D(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
}

inline const Vector3D PlaneSurface::binningPosition(BinningValue) const
{
  return center();
}

inline double
PlaneSurface::pathCorrection(const Vector3D&, const Vector3D& mom) const
{
  /// we can ignore the global position here
  return 1. / std::abs(normal().dot(mom.unit()));
}

inline Intersection
PlaneSurface::intersectionEstimate(const Vector3D&      gpos,
                                   const Vector3D&      gmom,
                                   NavigationDirection  navDir,
                                   const BoundaryCheck& bcheck,
                                   CorrFnc              correct) const
{
  // minimize the call to transform()
  const auto& tMatrix = transform().matrix();
  Vector3D    pnormal = tMatrix.block<3, 1>(0, 2).transpose();
  Vector3D    pcenter = tMatrix.block<3, 1>(0, 3).transpose();
  // position and direciton information
  Vector3D lpos = gpos;
  Vector3D ldir = gmom.unit();
  Vector3D solution(0., 0., 0.);
  double   path = std::numeric_limits<double>::infinity();
  // lemma : the solver -> should catch current values
  auto solve =
      [&solution, &path, &lpos, &ldir, &pnormal, &pcenter, &navDir]() -> bool {
    double denom = ldir.dot(pnormal);
    if (denom) {
      path     = (pnormal.dot((pcenter - lpos))) / (denom);
      solution = (lpos + path * ldir);
    }
    // is valid if it goes into the right direction
    return (!navDir || path * navDir >= 0.);
  };
  // solve first
  bool valid = solve();
  // if configured to correct, do it and solve again
  if (valid and correct and correct(lpos, ldir, path)) valid = solve();
  // evaluate (if necessary in terms of boundaries)
  // @todo: speed up isOnSurface - we know that it is on surface
  //  all we need is to check if it's inside bounds
  valid = bcheck ? (valid && isOnSurface(solution, bcheck)) : valid;
  // return the result
  return Intersection(solution, path, valid);
}

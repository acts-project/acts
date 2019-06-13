// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Vertexing/VertexingError.hpp"


template <typename bfield_t, typename input_track_t, typename propagator_t>
double Acts::ImpactPoint3dEstimator<bfield_t, input_track_t, propagator_t>::
    calculateDistance(const BoundParameters& trkParams,
                      const Vector3D& refPos) const {

  // Normalized momentum
  const Vector3D normMomentum = trkParams.momentum().normalized();

  // Calculate the path length
  double pathLength = (refPos - trkParams.position()).dot(normMomentum);

  // Point of closest approach in 3D
  Vector3D closestApp = trkParams.position() + pathLength * normMomentum;
  // Difference vector
  Vector3D deltaR = refPos - closestApp;

  // Return distance
  return deltaR.norm();
}


template <typename bfield_t, typename input_track_t, typename propagator_t>
Acts::Result<std::unique_ptr<const Acts::BoundParameters>>
Acts::ImpactPoint3dEstimator<bfield_t, input_track_t, propagator_t>::
    getParamsAtIP3d(const GeometryContext& gctx,
                    const MagneticFieldContext& mctx,
                    const BoundParameters& trkParams,
                    const Vector3D& vtxPos) const {
  Vector3D trkSurfaceCenter = trkParams.referenceSurface().center(gctx);

  double d0 = trkParams.parameters()[ParID_t::eLOC_D0];
  double z0 = trkParams.parameters()[ParID_t::eLOC_Z0];
  double phi = trkParams.parameters()[ParID_t::ePHI];
  double theta = trkParams.parameters()[ParID_t::eTHETA];
  double qOvP = trkParams.parameters()[ParID_t::eQOP];

  double cotTheta = 1. / std::tan(theta);
  double bZ = m_cfg.bField.getField(trkSurfaceCenter)[eZ];

  // The radius
  double r;
  // Curvature is infinite w/o b field
  if (bZ == 0. || std::abs(qOvP) < 1.e-15) {
    r = std::numeric_limits<double>::max();
  } else {
    // signed(!) r
    r = std::sin(theta) * units::Nat2SI<units::MOMENTUM>(1 / qOvP) / bZ;
  }

  Vector3D vec0 = trkSurfaceCenter + Vector3D(-(d0 - r) * std::sin(phi),
                                              (d0 - r) * std::cos(phi),
                                              z0 + r * phi * cotTheta);

  // Perform newton approximation method
  // this will change the value of phi
  auto res = performNewtonApproximation(vec0, vtxPos, phi, theta, r);
  if (!res.ok()) {
    return res.error();
  }

  Vector3D momDir(std::sin(theta) * std::cos(phi),
                  std::sin(phi) * std::sin(theta), std::cos(theta));

  // point of closest approach in 3D
  Vector3D pointCA3d =
      vec0 + r * Vector3D(-std::sin(phi), std::cos(phi), -cotTheta * phi);
  Vector3D deltaR = (pointCA3d - vtxPos).normalized();

  // corrected deltaR for small deviations from orthogonality
  Vector3D corrDeltaR = deltaR - (deltaR.dot(momDir)) * momDir;

  // get perpendicular direction vector
  Vector3D perpDir = momDir.cross(corrDeltaR);

  Transform3D thePlane;
  // rotation matrix
  thePlane.matrix().block(0, 0, 3, 1) = corrDeltaR;
  thePlane.matrix().block(0, 1, 3, 1) = perpDir;
  thePlane.matrix().block(0, 2, 3, 1) = momDir;
  // translation
  thePlane.matrix().block(0, 3, 3, 1) = vtxPos;

  std::shared_ptr<PlaneSurface> planeSurface =
      Surface::makeShared<PlaneSurface>(
          std::make_shared<Transform3D>(thePlane));

  PropagatorOptions pOptions(gctx, mctx);

  // Do the propagation to linPointPos
  auto result = m_cfg.propagator.propagate(trkParams, *planeSurface, pOptions);
  if (result.ok()) {
    return std::move((*result).endParameters);
  } else {
    return result.error();
  }
}

template <typename bfield_t, typename input_track_t, typename propagator_t>
Acts::Result<void> Acts::ImpactPoint3dEstimator<
    bfield_t, input_track_t,
    propagator_t>::performNewtonApproximation(const Vector3D& trkPos,
                                              const Vector3D& vtxPos,
                                              double& newPhi, double theta,
                                              double r) const {
  double sinNewPhi = -std::sin(newPhi);
  double cosNewPhi = std::cos(newPhi);

  int nIter = 0;
  bool hasConverged = false;

  double cotTheta = 1. / std::tan(theta);

  // start iteration until convergence or max iteration reached
  while (!hasConverged && nIter < m_cfg.maxIterations) {
    double x0 = trkPos.x();
    double y0 = trkPos.y();
    double z0 = trkPos.z();

    double xc = vtxPos.x();
    double yc = vtxPos.y();
    double zc = vtxPos.z();

    double derivative = (x0 - xc) * (-r * cosNewPhi) +
                        (y0 - yc) * r * sinNewPhi +
                        (z0 - zc - r * newPhi * cotTheta) * (-r * cotTheta);
    double secDerivative = r * (-(x0 - xc) * sinNewPhi - (y0 - yc) * cosNewPhi +
                                r * cotTheta * cotTheta);

    if (secDerivative < 0.) {
      return VertexingError::NumericFailure;
    }

    double deltaPhi = -derivative / secDerivative;

    newPhi += deltaPhi;
    sinNewPhi = -std::sin(newPhi);
    cosNewPhi = std::cos(newPhi);

    nIter += 1;

    if (std::abs(deltaPhi) < m_cfg.precision) {
      hasConverged = true;
    }
  }  // end while loop

  if (!hasConverged) {
    // max iterations reached but did not converge
    return VertexingError::NotConverged;
  }
  return {};
}


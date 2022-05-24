// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/SpacePointUtility.hpp"

#include <iostream>
namespace Acts {

std::pair<Vector2, SymMatrix2> SpacePointUtility::getLocalPosCov(
    const Measurement& meas) const {
  return std::visit(
      [](const auto& x) {
        auto expander = x.expander();
        BoundVector par = expander * x.parameters();
        BoundSymMatrix cov = expander * x.covariance() * expander.transpose();
        Vector2 lpar(par[BoundIndices::eBoundLoc0],
                     par[BoundIndices::eBoundLoc1]);
        SymMatrix2 lcov = cov.block<2, 2>(eBoundLoc0, eBoundLoc0);
        return std::make_pair(lpar, lcov);
      },
      meas);
}

std::pair<Vector3, Vector2> SpacePointUtility::globalCoords(
    const GeometryContext& gctx, const Measurement& meas) const {
  const auto& slink =
      std::visit([](const auto& x) { return &x.sourceLink(); }, meas);

  const auto geoId = slink->geometryId();

  const Surface* surface = m_config.trackingGeometry->findSurface(geoId);
  auto [localPos, localCov] = std::visit(
      [](const auto& measurement) {
        auto expander = measurement.expander();
        BoundVector par = expander * measurement.parameters();
        BoundSymMatrix cov =
            expander * measurement.covariance() * expander.transpose();
        // extract local position
        Vector2 lpar(par[eBoundLoc0], par[eBoundLoc1]);
        // extract local position covariance.
        SymMatrix2 lcov = cov.block<2, 2>(eBoundLoc0, eBoundLoc0);
        return std::make_pair(lpar, lcov);
      },
      meas);
  Vector3 globalPos = surface->localToGlobal(gctx, localPos, Vector3());
  RotationMatrix3 rotLocalToGlobal =
      surface->referenceFrame(gctx, globalPos, Vector3());

  // the space point requires only the variance of the transverse and
  // longitudinal position. reduce computations by transforming the
  // covariance directly from local to rho/z.
  //
  // compute Jacobian from global coordinates to rho/z
  //
  //         rho = sqrt(x² + y²)
  // drho/d{x,y} = (1 / sqrt(x² + y²)) * 2 * {x,y}
  //             = 2 * {x,y} / r
  //       dz/dz = 1 (duuh!)
  //
  auto x = globalPos[ePos0];
  auto y = globalPos[ePos1];
  auto scale = 2 / std::hypot(x, y);
  ActsMatrix<2, 3> jacXyzToRhoZ = ActsMatrix<2, 3>::Zero();
  jacXyzToRhoZ(0, ePos0) = scale * x;
  jacXyzToRhoZ(0, ePos1) = scale * y;
  jacXyzToRhoZ(1, ePos2) = 1;
  // compute Jacobian from local coordinates to rho/z
  ActsMatrix<2, 2> jac =
      jacXyzToRhoZ * rotLocalToGlobal.block<3, 2>(ePos0, ePos0);
  // compute rho/z variance
  ActsVector<2> var = (jac * localCov * jac.transpose()).diagonal();

  auto gcov = Vector2(var[0], var[1]);
  return std::make_pair(globalPos, gcov);
}

Vector2 SpacePointUtility::calcGlobalVars(const GeometryContext& gctx,
                                          const Measurement& measFront,
                                          const Measurement& measBack,
                                          const double theta) const {
  const auto var1 = getLoc0Var(measFront);
  const auto var2 = getLoc0Var(measBack);
  // strip1 and strip2 are tilted at +/- theta/2

  double sigma_x = std::hypot(var1, var2) / (2 * sin(theta * 0.5));
  double sigma_y = std::hypot(var1, var2) / (2 * cos(theta * 0.5));

  // projection to the surface with strip1.
  double sig_x1 = sigma_x * cos(0.5 * theta) + sigma_y * sin(0.5 * theta);
  double sig_y1 = sigma_y * cos(0.5 * theta) + sigma_x * sin(0.5 * theta);
  SymMatrix2 lcov;
  lcov << sig_x1, 0, 0, sig_y1;

  auto [localPos, localCov] = getLocalPosCov(measFront);

  const auto& slink_meas1 =
      std::visit([](const auto& x) { return &x.sourceLink(); }, measFront);

  const auto geoId = slink_meas1->geometryId();

  auto gcov = globalCov(gctx, geoId, localPos, lcov);

  return gcov;
}

double SpacePointUtility::getLoc0Var(const Measurement& meas) const {
  auto cov = std::visit(
      [](const auto& x) {
        auto expander = x.expander();
        BoundSymMatrix bcov = expander * x.covariance() * expander.transpose();
        SymMatrix2 lcov = bcov.block<2, 2>(eBoundLoc0, eBoundLoc0);
        return lcov;
      },
      meas);
  return cov(0, 0);
}

Vector2 SpacePointUtility::globalCov(const GeometryContext& gctx,
                                     const GeometryIdentifier& geoId,
                                     const Vector2& localPos,
                                     const SymMatrix2& localCov) const {
  Vector3 globalFakeMom(1, 1, 1);

  const Surface* surface = m_config.trackingGeometry->findSurface(geoId);

  Vector3 globalPos = surface->localToGlobal(gctx, localPos, globalFakeMom);
  RotationMatrix3 rotLocalToGlobal =
      surface->referenceFrame(gctx, globalPos, globalFakeMom);

  auto x = globalPos[ePos0];
  auto y = globalPos[ePos1];
  auto scale = 2 / std::hypot(x, y);
  ActsMatrix<2, 3> jacXyzToRhoZ = ActsMatrix<2, 3>::Zero();
  jacXyzToRhoZ(0, ePos0) = scale * x;
  jacXyzToRhoZ(0, ePos1) = scale * y;
  jacXyzToRhoZ(1, ePos2) = 1;
  // compute Jacobian from local coordinates to rho/z
  ActsMatrix<2, 2> jac =
      jacXyzToRhoZ * rotLocalToGlobal.block<3, 2>(ePos0, ePos0);
  // compute rho/z variance
  ActsVector<2> var = (jac * localCov * jac.transpose()).diagonal();

  auto gcov = Vector2(var[0], var[1]);

  return gcov;
}

const SourceLink* SpacePointUtility::getSourceLink(
    const Measurement meas) const {
  const SourceLink* slink =
      std::visit([](const auto& x) { return &x.sourceLink(); }, meas);
  return slink;
}

}  // namespace Acts

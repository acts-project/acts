// This file is part of the Acts project.
//
// Copyright (C) 2018-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename spacepoint_t, typename source_link_t>
Acts::Vector2
Acts::SingleHitSpacePointBuilder<spacepoint_t, source_link_t>::localCoords(
    const BoundVariantMeasurement<source_link_t>& meas) const {
  // Local position information
  auto par = meas.parameters();
  Acts::Vector2 local(par[Acts::BoundIndices::eBoundLoc0],
                      par[Acts::BoundIndices::eBoundLoc1]);
  return local;
}

template <typename spacepoint_t, typename source_link_t>
Acts::SingleHitSpacePointBuilder<spacepoint_t, source_link_t>::
    SingleHitSpacePointBuilder(Acts::SpacePointBuilderConfig cfg)
    : m_config(cfg) {}

template <typename spacepoint_t, typename source_link_t>
std::pair<Acts::Vector3, Acts::Vector2>
Acts::SingleHitSpacePointBuilder<spacepoint_t, source_link_t>::globalCoords(
    const GeometryContext& gctx,
    const BoundVariantMeasurement<source_link_t>& meas) const {
  auto slink = std::visit([](const auto& x) { return x.sourceLink(); }, meas);

  const auto geoId = slink.geometryId();
  const Acts::Surface* surface = m_config.trackingGeometry->findSurface(geoId);

  auto [localPos, localCov] = std::visit(
      [](const auto& meas) {
        auto expander = meas.expander();
        Acts::BoundVector par = expander * meas.parameters();
        Acts::BoundSymMatrix cov =
            expander * meas.covariance() * expander.transpose();
        // extract local position
        Acts::Vector2 lpar(par[Acts::eBoundLoc0], par[Acts::eBoundLoc1]);
        // extract local position covariance.
        Acts::SymMatrix2 lcov =
            cov.block<2, 2>(Acts::eBoundLoc0, Acts::eBoundLoc0);
        return std::make_pair(lpar, lcov);
      },
      meas);

  // transform local position to global coordinates
  Acts::Vector3 globalFakeMom(1, 1, 1);
  Acts::Vector3 globalPos =
      surface->localToGlobal(gctx, localPos, globalFakeMom);
  Acts::RotationMatrix3 rotLocalToGlobal =
      surface->referenceFrame(gctx, globalPos, globalFakeMom);

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
  auto x = globalPos[Acts::ePos0];
  auto y = globalPos[Acts::ePos1];
  auto scale = 2 / std::hypot(x, y);
  Acts::ActsMatrix<2, 3> jacXyzToRhoZ = Acts::ActsMatrix<2, 3>::Zero();
  jacXyzToRhoZ(0, Acts::ePos0) = scale * x;
  jacXyzToRhoZ(0, Acts::ePos1) = scale * y;
  jacXyzToRhoZ(1, Acts::ePos2) = 1;
  // compute Jacobian from local coordinates to rho/z
  Acts::ActsMatrix<2, 2> jac =
      jacXyzToRhoZ * rotLocalToGlobal.block<3, 2>(Acts::ePos0, Acts::ePos0);
  // compute rho/z variance
  Acts::ActsVector<2> var = (jac * localCov * jac.transpose()).diagonal();

  auto gcov = Acts::Vector2(var[0], var[1]);
  return std::make_pair(globalPos, gcov);
}

template <typename spacepoint_t, typename source_link_t>
void Acts::SingleHitSpacePointBuilder<spacepoint_t, source_link_t>::
    calculateSpacePoints(
        const GeometryContext& gctx,
        const std::vector<BoundVariantMeasurement<source_link_t>>& measurements,
        std::vector<spacepoint_t>& spacePointStorage) const {
  // Set the space point for all stored hits
  for (const auto& meas : measurements) {
    auto [gPos, gCov] = globalCoords(gctx, meas);
    auto slink = std::visit([](const auto& x) { return x.sourceLink(); }, meas);

    spacePointStorage.emplace_back(gPos, gCov[0], gCov[1], slink.index());
  }
}
// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

ActsExamples::SpacePointMaker::SpacePointMaker(Config cfg,
                                               Acts::Logging::Level lvl)
    : BareAlgorithm("SpacePointMaker", lvl), m_cfg(std::move(cfg)) {
  if (m_cfg.inputSourceLinks.empty()) {
    throw std::invalid_argument("Missing source link input collection");
  }
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurement input collection");
  }
  if (m_cfg.outputSpacePoints.empty()) {
    throw std::invalid_argument("Missing space point output collection");
  }
  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }
  if (m_cfg.geometrySelection.empty()) {
    throw std::invalid_argument("Missing space point maker geometry selection");
  }
  // ensure geometry selection contains only valid inputs
  for (const auto& geoId : m_cfg.geometrySelection) {
    if ((geoId.approach() != 0u) or (geoId.boundary() != 0u) or
        (geoId.sensitive() != 0u)) {
      throw std::invalid_argument(
          "Invalid geometry selection: only volume and layer are allowed to be "
          "set");
    }
  }
  // remove geometry selection duplicates
  //
  // the geometry selections must be mutually exclusive, i.e. if we have a
  // selection that contains both a volume and a layer within that same volume,
  // we would create the space points for the layer twice.
  auto isDuplicate = [](Acts::GeometryIdentifier ref,
                        Acts::GeometryIdentifier cmp) {
    // code assumes ref < cmp and that only volume and layer can be non-zero
    // root node always contains everything
    if (ref.volume() == 0) {
      return true;
    }
    // unequal volumes always means separate hierarchies
    if (ref.volume() != cmp.volume()) {
      return false;
    }
    // within the same volume hierarchy only consider layers
    return (ref.layer() == cmp.layer());
  };
  auto geoSelBeg = m_cfg.geometrySelection.begin();
  auto geoSelEnd = m_cfg.geometrySelection.end();
  // sort geometry selection so the unique filtering works
  std::sort(geoSelBeg, geoSelEnd);
  auto geoSelLastUnique = std::unique(geoSelBeg, geoSelEnd, isDuplicate);
  if (geoSelLastUnique != geoSelEnd) {
    ACTS_WARNING("Removed " << std::distance(geoSelLastUnique, geoSelEnd)
                            << " geometry selection duplicates");
    m_cfg.geometrySelection.erase(geoSelLastUnique, geoSelEnd);
  }
  ACTS_INFO("Space point geometry selection:");
  for (const auto& geoId : m_cfg.geometrySelection) {
    ACTS_INFO("  " << geoId);
  }
}

ActsExamples::ProcessCode ActsExamples::SpacePointMaker::execute(
    const AlgorithmContext& ctx) const {
  const auto& sourceLinks =
      ctx.eventStore.get<IndexSourceLinkContainer>(m_cfg.inputSourceLinks);
  const auto& measurements =
      ctx.eventStore.get<MeasurementContainer>(m_cfg.inputMeasurements);

  SimSpacePointContainer spacePoints;
  spacePoints.reserve(sourceLinks.size());

  for (Acts::GeometryIdentifier geoId : m_cfg.geometrySelection) {
    // select volume/layer depending on what is set in the geometry id
    auto range = selectLowestNonZeroGeometryObject(sourceLinks, geoId);
    // groupByModule only works with geometry containers, not with an
    // arbitrary range. do the equivalent grouping manually
    auto groupedByModule = makeGroupBy(range, detail::GeometryIdGetter());

    for (auto [moduleGeoId, moduleSourceLinks] : groupedByModule) {
      // find corresponding surface
      const Acts::Surface* surface =
          m_cfg.trackingGeometry->findSurface(moduleGeoId);
      if (surface == nullptr) {
        ACTS_ERROR("Could not find surface " << moduleGeoId);
        return ProcessCode::ABORT;
      }

      for (auto& sourceLink : moduleSourceLinks) {
        // extract a local position/covariance independent from the concrecte
        // measurement content. since we do not know if and where the local
        // parameters are contained in the measurement parameters vector, they
        // are transformed to the bound space where we do know their location.
        // if the local parameters are not measured, this results in a
        // zero location, which is a reasonable default fall-back.
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
            measurements[sourceLink.get().index()]);

        // transform local position to global coordinates
        Acts::Vector3 globalFakeMom(1, 1, 1);
        Acts::Vector3 globalPos =
            surface->localToGlobal(ctx.geoContext, localPos, globalFakeMom);
        Acts::RotationMatrix3 rotLocalToGlobal =
            surface->referenceFrame(ctx.geoContext, globalPos, globalFakeMom);

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
            jacXyzToRhoZ *
            rotLocalToGlobal.block<3, 2>(Acts::ePos0, Acts::ePos0);
        // compute rho/z variance
        Acts::ActsVector<2> var = (jac * localCov * jac.transpose()).diagonal();

        // construct space point in global coordinates
        spacePoints.emplace_back(globalPos, var[0], var[1],
                                 sourceLink.get().index());
      }
    }
  }
  spacePoints.shrink_to_fit();

  ACTS_DEBUG("Created " << spacePoints.size() << " space points");
  ctx.eventStore.add(m_cfg.outputSpacePoints, std::move(spacePoints));

  return ActsExamples::ProcessCode::SUCCESS;
}

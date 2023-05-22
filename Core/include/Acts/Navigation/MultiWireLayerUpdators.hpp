// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Detector/detail/GridAxisGenerators.hpp"
#include "Acts/Detector/detail/IndexedGridFiller.hpp"
#include "Acts/Detector/detail/IndexedSurfacesGenerator.hpp"
#include "Acts/Detector/detail/ReferenceGenerators.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Navigation/NavigationStateFillers.hpp"
#include "Acts/Navigation/NavigationStateUpdators.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

#include <memory>
#include <tuple>

namespace Acts {
namespace Experimental {

using grid_type = detail::GridAxisGenerators::EqBoundEqBound::grid_type<
    std::vector<std::size_t>>;
using IndexedSurfacesGridType = IndexedSurfacesImpl<grid_type>;

struct MultiWireLayerImpl : public INavigationDelegate {
  IndexedSurfacesGridType getIndexedSurfacesFromVolume(
      NavigationState& nState) const {
    auto gridSurfaces = nState.currentVolume->surfaces();
    GeometryContext gctx;

    std::array<std::pair<float, float>, 3> min_max_multilayer;
    std::fill(min_max_multilayer.begin(), min_max_multilayer.end(),
              std::make_pair<float, float>(std::numeric_limits<float>::max(),
                                           -std::numeric_limits<float>::max()));

    for (auto& gridSurf : gridSurfaces) {
      min_max_multilayer[0].first = std::min(min_max_multilayer[0].first,
                                             (float)gridSurf->center(gctx).x());
      min_max_multilayer[0].second = std::max(
          min_max_multilayer[0].second, (float)gridSurf->center(gctx).x());

      min_max_multilayer[1].first = std::min(min_max_multilayer[1].first,
                                             (float)gridSurf->center(gctx).y());
      min_max_multilayer[1].second = std::max(
          min_max_multilayer[1].second, (float)gridSurf->center(gctx).y());

      min_max_multilayer[2].first = std::min(min_max_multilayer[2].first,
                                             (float)gridSurf->center(gctx).z());
      min_max_multilayer[2].second = std::max(
          min_max_multilayer[2].second, (float)gridSurf->center(gctx).z());
    }
    auto radius = gridSurfaces.front()->bounds().values()[0];

    float binWidthX = 2 * radius, binWidthY = sqrt(3) * radius;

    float rangeX1 = min_max_multilayer[2].first - radius,
          rangeX2 = min_max_multilayer[2].second + radius;
    float rangeY1 = min_max_multilayer[1].first - radius,
          rangeY2 = min_max_multilayer[1].second + radius;

    std::size_t nBinsX =
                    static_cast<std::size_t>((rangeX2 - rangeX1) / binWidthX),
                nBinsY =
                    static_cast<std::size_t>((rangeY2 - rangeY1) / binWidthY);

    detail::IndexedSurfacesGenerator<decltype(gridSurfaces)> irSurfaces_grid{
        gridSurfaces, {}, {binZ, binY}, {1u, 0u}};

    detail::GridAxisGenerators::EqBoundEqBound aGenerator_grid{
        {rangeX1, rangeX2}, nBinsX, {rangeY1, rangeY2}, nBinsY};

    detail::CenterReferenceGenerator rGenerator;

    auto indexedSurfaces = irSurfaces_grid(gctx, aGenerator_grid, rGenerator);

    using GridType =
        decltype(aGenerator_grid)::grid_type<std::vector<std::size_t>>;
    using DelegateType = IndexedSurfacesAllPortalsImpl<GridType>;

    const auto* instance_multilayer = indexedSurfaces.instance();
    auto castedDelegate_grid =
        dynamic_cast<const DelegateType*>(instance_multilayer);

    const auto& chainedUpdators_grid = castedDelegate_grid->updators;
    const auto& indexedSurfaces_multilayer =
        std::get<IndexedSurfacesImpl<GridType>>(chainedUpdators_grid);

    return indexedSurfaces_multilayer;
  }

  inline void update(const Acts::GeometryContext& gctx,
                     Acts::Experimental::NavigationState& nState) const {
    if (nState.currentVolume == nullptr) {
      throw std::runtime_error(
          "MultiWireLayerImpl: no detector volume set to navigation state.");
    }

    // A volume switch has happened, update list of portals & surfaces

    if (nState.surfaceCandidates.empty()) {
      // Fill internal portals if existing
      for (const auto v : nState.currentVolume->volumes()) {
        const auto& iPortals = v->portals();
        Acts::Experimental::PortalsFiller::fill(nState, iPortals);
      }

      // Fill with the indexed surfaces and the portals of the volume

      const auto m_indexedSurfaces = getIndexedSurfacesFromVolume(nState);

      auto tposition = nState.position;
      auto tdirection = nState.direction;

      const auto& grid = m_indexedSurfaces.grid;
      const auto& extractor = m_indexedSurfaces.extractor;
      auto nBins = grid.numLocalBins();
      auto binWidth = grid.binWidth();
      auto minPos = grid.minPosition();
      auto maxPos = grid.maxPosition();

      Vector3 origin = Vector3(0, minPos[1], 0);
      Translation3 trans(origin);
      Transform3 transform(trans);

      auto dy = binWidth[1];
      auto phi = Acts::VectorHelpers::phi(tdirection);

      BinUtility binUtility2D(nBins[0], minPos[0], maxPos[0], open, binX,
                              transform);
      binUtility2D +=
          BinUtility(nBins[1], minPos[1], maxPos[1], open, binY, transform);

      for (std::size_t i = 0; i < nBins[1]; i++) {
        // auto lbin = grid.localBinsFromPosition(tposition);
        // auto lbin = grid.globalBinFromPosition(tposition);
        Vector2 localposition =
            Vector2(tposition.x() + tolerance, tposition.y() + tolerance);
        /* std::cout<<"local
         position="<<localposition.x()<<","<<localposition.y()<<std::endl;
         std::cout<<"bin utility bin="<<binUtility2D.bin(tposition)<<std::endl;
         std::cout<<"bin utility
         bin="<<binUtility2D.bin(tposition,1)<<std::endl;*/

        std::array<size_t, 2> lbin = {binUtility2D.bin(localposition) + 1,
                                      binUtility2D.bin(localposition, 1) + 1};
        auto extracted =
            extractor.extract(gctx, nState, grid.atLocalBins(lbin));
        // auto extracted = extractor.extract(gctx, nState,
        // grid.atPosition(tposition));
        std::cout << "lbin=" << lbin[0] << "," << lbin[1] << std::endl;
        std::cout << "extracted size=" << extracted.size() << std::endl;
        std::cout << "extracted centers:" << std::endl;
        std::cout << extracted[0]->center(gctx).y() << "\t";
        std::cout << extracted[1]->center(gctx).y() << "\t";
        std::cout << extracted[2]->center(gctx).y() << "\t";

        Acts::Experimental::SurfacesFiller::fill(nState, extracted);

        tposition.x() = tposition.x() + dy / tan(phi);
        tposition.y() = tposition.y() + dy;
        /*std::cout<<"in the MultiWireLayer"<<std::endl;
        std::cout<<nState.surfaceCandidates.size()<<std::endl;
        std::cout<<"is inside="<<binUtility2D.inside(localposition)<<std::endl;
        std::cin.ignore();*/
      }

      const auto& portals = nState.currentVolume->portals();
      Acts::Experimental::PortalsFiller::fill(nState, portals);
    }

    updateCandidates(gctx, nState);
  }

 private:
  float tolerance = 0.5;
};

inline static Acts::Experimental::SurfaceCandidatesUpdator MultiWireLayer() {
  auto ap = std::make_unique<const MultiWireLayerImpl>();
  SurfaceCandidatesUpdator nStateUpdator;
  nStateUpdator.connect<&MultiWireLayerImpl::update>(std::move(ap));
  return nStateUpdator;
}

}  // namespace Experimental
}  // namespace Acts

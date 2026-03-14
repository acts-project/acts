// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/SurfaceArray.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Ranges.hpp"

#include <format>
#include <ranges>
#include <utility>

using namespace Acts::Ranges;
namespace Acts {

// implementation for pure virtual destructor of ISurfaceGridLookup
SurfaceArray::ISurfaceGridLookup::~ISurfaceGridLookup() = default;

SurfaceArray::SurfaceArray(std::unique_ptr<ISurfaceGridLookup> gridLookup,
                           std::vector<std::shared_ptr<const Surface>> surfaces,
                           const Transform3& transform)
    : p_gridLookup(std::move(gridLookup)),
      m_surfaces(std::move(surfaces)),
      m_surfacesRawPointers(unpackSmartPointers(m_surfaces)),
      m_transform(transform) {
  if (p_gridLookup != nullptr) {
    if (const auto& grid = p_gridLookup->getGridView()) {
      checkGrid(grid.value());
    }
  }
}

SurfaceArray::SurfaceArray(std::shared_ptr<const Surface> srf)
    : p_gridLookup(std::make_unique<SingleElementLookup>(srf.get())),
      m_surfaces({std::move(srf)}) {
  m_surfacesRawPointers.push_back(m_surfaces.at(0).get());
}

std::ostream& SurfaceArray::toStream(const GeometryContext& /*gctx*/,
                                     std::ostream& sl) const {
  sl << std::fixed << std::setprecision(4);
  sl << "SurfaceArray:" << std::endl;
  sl << " - no surfaces: " << m_surfaces.size() << std::endl;

  auto axes = p_gridLookup->getAxes();

  for (std::size_t j = 0; j < axes.size(); ++j) {
    AxisBoundaryType bdt = axes.at(j)->getBoundaryType();
    sl << " - axis " << (j + 1) << std::endl;
    sl << "   - boundary type: ";
    if (bdt == AxisBoundaryType::Open) {
      sl << "open";
    }
    if (bdt == AxisBoundaryType::Bound) {
      sl << "bound";
    }
    if (bdt == AxisBoundaryType::Closed) {
      sl << "closed";
    }
    sl << std::endl;
    sl << "   - type: "
       << (axes.at(j)->isEquidistant() ? "equidistant" : "variable")
       << std::endl;
    sl << "   - n bins: " << axes.at(j)->getNBins() << std::endl;
    sl << "   - bin edges: [ ";
    auto binEdges = axes.at(j)->getBinEdges();
    for (std::size_t i = 0; i < binEdges.size(); ++i) {
      if (i > 0) {
        sl << ", ";
      }
      auto binEdge = binEdges.at(i);
      // Do not display negative zeroes
      sl << ((std::abs(binEdge) >= 5e-4) ? binEdge : 0.0);
    }
    sl << " ]" << std::endl;
  }
  return sl;
}

void SurfaceArray::checkGrid(AnyGridConstView<SurfaceVector> grid) {
  std::set allSurfaces =
      m_surfaces |
      std::views::transform([](const auto& sp) { return sp.get(); }) |
      to<std::set>;
  std::set<const Surface*> seenSurface;
  auto bins = grid.numLocalBins();
  for (std::size_t i = 0; i <= bins.at(0); ++i) {
    for (std::size_t j = 0; j <= bins.at(1); ++j) {
      const auto& surfaces = grid.atLocalBins({i, j});
      for (const auto& srf : surfaces) {
        seenSurface.insert(srf);
      }
    }
  }

  if (allSurfaces != seenSurface) {
    std::set<const Surface*> diff;
    std::ranges::set_difference(allSurfaces, seenSurface,
                                std::inserter(diff, diff.begin()));

    throw std::logic_error(
        std::format("SurfaceArray grid does not contain all surfaces provided! "
                    "{} surfaces not seen",
                    diff.size()));
  }
}

}  // namespace Acts

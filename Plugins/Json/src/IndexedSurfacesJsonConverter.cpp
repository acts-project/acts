// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/IndexedSurfacesJsonConverter.hpp"

#include "Acts/Detector/detail/GridAxisGenerators.hpp"
#include "Acts/Navigation/NavigationStateUpdators.hpp"
#include "Acts/Plugins/Json/GridJsonConverter.hpp"
#include "Acts/Plugins/Json/UtilitiesJsonConverter.hpp"

#include <array>
#include <memory>
#include <tuple>
#include <vector>

namespace {

/// @brief  Helper function to create and connect the IndexedSurfacesImpl
///
/// @tparam grid_type the type of the grid, indicates also the dimension
///
/// @param grid the grid object
/// @param bv the bin value array
/// @param transform the transform for the indexed surfaces inmplementaiton
///
/// @return a connected SurfaceCandidatesUpdator object
template <typename grid_type>
Acts::Experimental::SurfaceCandidatesUpdator createUpdator(
    grid_type&& grid, const std::array<Acts::BinningValue, grid_type::DIM>& bv,
    const Acts::Transform3& transform) {
  Acts::Experimental::IndexedSurfacesImpl<grid_type> indexedSurfaces(
      std::move(grid), bv, transform);

  // The portal delegate
  Acts::Experimental::AllPortalsImpl allPortals;

  // The chained delegate: indexed surfaces and all portals
  using DelegateType = Acts::Experimental::IndexedSurfacesAllPortalsImpl<
      grid_type, Acts::Experimental::IndexedSurfacesImpl>;
  auto indexedSurfacesAllPortals = std::make_unique<const DelegateType>(
      std::tie(allPortals, indexedSurfaces));

  // Create the delegate and connect it
  Acts::Experimental::SurfaceCandidatesUpdator nStateUpdator;
  nStateUpdator.connect<&DelegateType::update>(
      std::move(indexedSurfacesAllPortals));

  return nStateUpdator;
}

}  // namespace

Acts::Experimental::SurfaceCandidatesUpdator
Acts::IndexedSurfacesJsonConverter::fromJson(
    const nlohmann::json& jSurfaceNavigation) {
  // The return object
  Experimental::SurfaceCandidatesUpdator sfCandidates;

  // Helper extractor
  auto eqExtractor = [](const nlohmann::json& jAxis)
      -> std::tuple<std::array<ActsScalar, 2u>, std::size_t> {
    std::array<ActsScalar, 2u> range = jAxis["range"];
    std::size_t bins = jAxis["bins"];
    return std::make_tuple(range, bins);
  };

  auto vExtractor = [](const nlohmann::json& jAxis) -> std::vector<ActsScalar> {
    return std::vector<ActsScalar>(jAxis["boundaries"]);
  };

  // Peek into the json object to understand what to do
  if (jSurfaceNavigation["type"] == "IndexedSurfaces") {
    if (jSurfaceNavigation.find("grid") != jSurfaceNavigation.end()) {
      Transform3 transform =
          Transform3JsonConverter::fromJson(jSurfaceNavigation["transform"]);
      auto jGrid = jSurfaceNavigation["grid"];
      auto jCasts =
          jSurfaceNavigation["casts"].get<std::vector<BinningValue>>();
      auto jAxes = jGrid["axes"];

      // 1D cases
      if (jAxes.size() == 1u) {
        BinningValue bValue = jCasts[0u];
        auto jAxis = jAxes[0u];

        detail::AxisType axisType = jAxis["type"];
        detail::AxisBoundaryType axisBoundaryType = jAxis["boundary_type"];

        // Equidistant axis
        if (axisType == detail::AxisType::Equidistant) {
          auto [range, bins] = eqExtractor(jAxis);
          if (axisBoundaryType == detail::AxisBoundaryType::Closed) {
            EqClosed ecAG{range, bins};
            auto grid = GridJsonConverter::fromJson(jGrid, ecAG);
            return createUpdator(std::move(grid), {bValue}, transform);
          } else {
            EqBound ebAG{range, bins};
            auto grid = GridJsonConverter::fromJson(jGrid, ebAG);
            return createUpdator(std::move(grid), {bValue}, transform);
          }
        } else {
          // Variable type
          if (axisBoundaryType == detail::AxisBoundaryType::Closed) {
            VarClosed vcAG{vExtractor(jAxis)};
            auto grid = GridJsonConverter::fromJson(jGrid, vcAG);
            return createUpdator(std::move(grid), {bValue}, transform);
          } else {
            VarBound vbAG{vExtractor(jAxis)};
            auto grid = GridJsonConverter::fromJson(jGrid, vbAG);
            return createUpdator(std::move(grid), {bValue}, transform);
          }
        }
      } else if (jAxes.size() == 2u) {
        // This currently writes out only the main options of 2D grids
        // nota bene: it assumes if one axis is closed, it is axis B

        BinningValue bValueA = jCasts[0u];
        BinningValue bValueB = jCasts[1u];
        auto jAxisA = jAxes[0u];
        auto jAxisB = jAxes[1u];

        detail::AxisType axisTypeA = jAxisA["type"];
        detail::AxisType axisTypeB = jAxisB["type"];
        detail::AxisBoundaryType axisBoundaryTypeB = jAxisB["boundary_type"];

        if (axisBoundaryTypeB != detail::AxisBoundaryType::Closed) {
          // First axis equidistant
          if (axisTypeA == detail::AxisType::Equidistant) {
            auto [rangeA, binsA] = eqExtractor(jAxisA);
            if (axisTypeB == detail::AxisType::Equidistant) {
              auto [rangeB, binsB] = eqExtractor(jAxisB);
              EqBoundEqBound ebebAG{rangeA, binsA, rangeB, binsB};
              auto grid = GridJsonConverter::fromJson(jGrid, ebebAG);
              return createUpdator(std::move(grid), {bValueA, bValueB},
                                   transform);
            } else {
              EqBoundVarBound ebvbAG{rangeA, binsA, vExtractor(jAxisB)};
              auto grid = GridJsonConverter::fromJson(jGrid, ebvbAG);
              return createUpdator(std::move(grid), {bValueA, bValueB},
                                   transform);
            }
          } else {
            if (axisTypeB == detail::AxisType::Equidistant) {
              auto [rangeB, binsB] = eqExtractor(jAxisB);
              VarBoundEqBound vbebAG{vExtractor(jAxisA), rangeB, binsB};
              auto grid = GridJsonConverter::fromJson(jGrid, vbebAG);
              return createUpdator(std::move(grid), {bValueA, bValueB},
                                   transform);
            } else {
              VarBoundVarBound vbvbAG{vExtractor(jAxisA), vExtractor(jAxisB)};
              auto grid = GridJsonConverter::fromJson(jGrid, vbvbAG);
              return createUpdator(std::move(grid), {bValueA, bValueB},
                                   transform);
            }
          }
        } else {
          // Closed cases
          if (axisTypeA == detail::AxisType::Equidistant) {
            auto [rangeA, binsA] = eqExtractor(jAxisA);
            if (axisTypeB == detail::AxisType::Equidistant) {
              auto [rangeB, binsB] = eqExtractor(jAxisB);
              EqBoundEqClosed ebecAG{rangeA, binsA, rangeB, binsB};
              auto grid = GridJsonConverter::fromJson(jGrid, ebecAG);
              return createUpdator(std::move(grid), {bValueA, bValueB},
                                   transform);
            } else {
              EqBoundVarClosed ebvcAG{rangeA, binsA, vExtractor(jAxisB)};
              auto grid = GridJsonConverter::fromJson(jGrid, ebvcAG);
              return createUpdator(std::move(grid), {bValueA, bValueB},
                                   transform);
            }
          } else {
            if (axisTypeB == detail::AxisType::Equidistant) {
              auto [rangeB, binsB] = eqExtractor(jAxisB);
              VarBoundEqClosed vbecAG{vExtractor(jAxisA), rangeB, binsB};
              auto grid = GridJsonConverter::fromJson(jGrid, vbecAG);
              return createUpdator(std::move(grid), {bValueA, bValueB},
                                   transform);
            } else {
              VarBoundVarClosed vbvcAG{vExtractor(jAxisA), vExtractor(jAxisB)};
              auto grid = GridJsonConverter::fromJson(jGrid, vbvcAG);
              return createUpdator(std::move(grid), {bValueA, bValueB},
                                   transform);
            }
          }
        }
      }
    }
  }
  // Return the object
  return Experimental::tryAllPortalsAndSurfaces();
}

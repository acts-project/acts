// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/detail/GridAxisGenerators.hpp"
#include "Acts/Plugins/Json/AlgebraJsonConverter.hpp"
#include "Acts/Plugins/Json/GridJsonConverter.hpp"
#include "Acts/Utilities/detail/AxisFwd.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

#include <tuple>

namespace Acts {
using namespace Experimental::detail::GridAxisGenerators;

namespace IndexedGridJsonHelper {

// Generate the possible axes in this case
static auto s_possibleAxes =
    std::tuple<EqBound, EqOpen, EqClosed,
               // All 1D Var  options
               VarBound, VarOpen, VarClosed,
               // All 2D EqEq options
               EqBoundEqBound, EqBoundEqOpen, EqBoundEqClosed, EqOpenEqBound,
               EqOpenEqOpen, EqOpenEqClosed, EqClosedEqBound, EqClosedEqOpen,
               EqClosedEqClosed,
               // All 2D EqVar options
               EqBoundVarBound, EqBoundVarOpen, EqBoundVarClosed,
               EqOpenVarBound, EqOpenVarOpen, EqOpenVarClosed, EqClosedVarBound,
               EqClosedVarOpen, EqClosedVarClosed,
               // All 2D VarEq options
               VarBoundEqBound, VarBoundEqOpen, VarBoundEqClosed,
               VarOpenEqBound, VarOpenEqOpen, VarOpenEqClosed, VarClosedEqBound,
               VarClosedEqOpen, VarClosedEqClosed,
               // All 2D VarEq options
               VarBoundVarBound, VarBoundVarOpen, VarBoundVarClosed,
               VarOpenVarBound, VarOpenVarOpen, VarOpenVarClosed,
               VarClosedVarBound, VarClosedVarOpen, VarClosedVarClosed>{};

/// @brief The actual conversion method
template <typename index_grid>
nlohmann::json convertImpl(const index_grid& indexGrid, bool detray = false) {
  nlohmann::json jIndexedGrid;

  // Fill the casts
  nlohmann::json jCasts;
  // 1D casts
  if constexpr (index_grid::grid_type::DIM == 1u) {
    jCasts.push_back(indexGrid.casts[0u]);
  }
  // 1D casts
  if constexpr (index_grid::grid_type::DIM == 2u) {
    jCasts.push_back(indexGrid.casts[0u]);
    jCasts.push_back(indexGrid.casts[1u]);
  }
  jIndexedGrid["casts"] = jCasts;
  jIndexedGrid["transform"] =
      Transform3JsonConverter::toJson(indexGrid.transform);
  if (detray) {
    jIndexedGrid = GridJsonConverter::toJsonDetray(indexGrid.grid);
  } else {
    jIndexedGrid["grid"] = GridJsonConverter::toJson(indexGrid.grid);
  }
  return jIndexedGrid;
}

/// @brief Helper method to be used for the creation of surface updators
/// and volume updates
///
/// @tparam updator_type
/// @tparam generator_type
///
/// @param jUpdator The corresponding json object
/// @param jIndicator the string indicator which one it is
///
/// @return the updator type
template <typename updator_type, typename generator_type>
updator_type generateFromJson(const nlohmann::json& jUpdator,
                              const std::string& jIndicator) {
  generator_type generator;

  using ValueType = typename generator_type::value_type;

  /// Helper extractor for equidistant axis
  /// @param jAxis is the axis
  auto eqExtractor = [](const nlohmann::json& jAxis)
      -> std::tuple<std::array<ActsScalar, 2u>, std::size_t> {
    std::array<ActsScalar, 2u> range = jAxis["range"];
    std::size_t bins = jAxis["bins"];
    return std::make_tuple(range, bins);
  };

  /// Helper extractor for variable axis
  /// @param jAxis the axis
  auto vExtractor = [](const nlohmann::json& jAxis) -> std::vector<ActsScalar> {
    return std::vector<ActsScalar>(jAxis["boundaries"]);
  };

  // Peek into the json object to understand what to do
  if (jUpdator["type"] == jIndicator) {
    if (jUpdator.find("grid") != jUpdator.end()) {
      Transform3 transform =
          Transform3JsonConverter::fromJson(jUpdator["transform"]);
      auto jGrid = jUpdator["grid"];
      auto jCasts = jUpdator["casts"].get<std::vector<BinningValue>>();
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
            auto grid =
                GridJsonConverter::fromJson<EqClosed, ValueType>(jGrid, ecAG);
            return generator.createUpdator(std::move(grid), {bValue},
                                           transform);
          } else {
            EqBound ebAG{range, bins};
            auto grid =
                GridJsonConverter::fromJson<EqBound, ValueType>(jGrid, ebAG);
            return generator.createUpdator(std::move(grid), {bValue},
                                           transform);
          }
        } else {
          // Variable type
          if (axisBoundaryType == detail::AxisBoundaryType::Closed) {
            VarClosed vcAG{vExtractor(jAxis)};
            auto grid =
                GridJsonConverter::fromJson<VarClosed, ValueType>(jGrid, vcAG);
            return generator.createUpdator(std::move(grid), {bValue},
                                           transform);
          } else {
            VarBound vbAG{vExtractor(jAxis)};
            auto grid =
                GridJsonConverter::fromJson<VarBound, ValueType>(jGrid, vbAG);
            return generator.createUpdator(std::move(grid), {bValue},
                                           transform);
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
              auto grid =
                  GridJsonConverter::fromJson<EqBoundEqBound, ValueType>(
                      jGrid, ebebAG);
              return generator.createUpdator(std::move(grid),
                                             {bValueA, bValueB}, transform);
            } else {
              EqBoundVarBound ebvbAG{rangeA, binsA, vExtractor(jAxisB)};
              auto grid =
                  GridJsonConverter::fromJson<EqBoundVarBound, ValueType>(
                      jGrid, ebvbAG);
              return generator.createUpdator(std::move(grid),
                                             {bValueA, bValueB}, transform);
            }
          } else {
            if (axisTypeB == detail::AxisType::Equidistant) {
              auto [rangeB, binsB] = eqExtractor(jAxisB);
              VarBoundEqBound vbebAG{vExtractor(jAxisA), rangeB, binsB};
              auto grid =
                  GridJsonConverter::fromJson<VarBoundEqBound, ValueType>(
                      jGrid, vbebAG);
              return generator.createUpdator(std::move(grid),
                                             {bValueA, bValueB}, transform);
            } else {
              VarBoundVarBound vbvbAG{vExtractor(jAxisA), vExtractor(jAxisB)};
              auto grid =
                  GridJsonConverter::fromJson<VarBoundVarBound, ValueType>(
                      jGrid, vbvbAG);
              return generator.createUpdator(std::move(grid),
                                             {bValueA, bValueB}, transform);
            }
          }
        } else {
          // Closed cases
          if (axisTypeA == detail::AxisType::Equidistant) {
            auto [rangeA, binsA] = eqExtractor(jAxisA);
            if (axisTypeB == detail::AxisType::Equidistant) {
              auto [rangeB, binsB] = eqExtractor(jAxisB);
              EqBoundEqClosed ebecAG{rangeA, binsA, rangeB, binsB};
              auto grid =
                  GridJsonConverter::fromJson<EqBoundEqClosed, ValueType>(
                      jGrid, ebecAG);
              return generator.createUpdator(std::move(grid),
                                             {bValueA, bValueB}, transform);
            } else {
              EqBoundVarClosed ebvcAG{rangeA, binsA, vExtractor(jAxisB)};
              auto grid =
                  GridJsonConverter::fromJson<EqBoundVarClosed, ValueType>(
                      jGrid, ebvcAG);
              return generator.createUpdator(std::move(grid),
                                             {bValueA, bValueB}, transform);
            }
          } else {
            if (axisTypeB == detail::AxisType::Equidistant) {
              auto [rangeB, binsB] = eqExtractor(jAxisB);
              VarBoundEqClosed vbecAG{vExtractor(jAxisA), rangeB, binsB};
              auto grid =
                  GridJsonConverter::fromJson<VarBoundEqClosed, ValueType>(
                      jGrid, vbecAG);
              return generator.createUpdator(std::move(grid),
                                             {bValueA, bValueB}, transform);
            } else {
              VarBoundVarClosed vbvcAG{vExtractor(jAxisA), vExtractor(jAxisB)};
              auto grid =
                  GridJsonConverter::fromJson<VarBoundVarClosed, ValueType>(
                      jGrid, vbvcAG);
              return generator.createUpdator(std::move(grid),
                                             {bValueA, bValueB}, transform);
            }
          }
        }
      }
    }
  }
  // The return object
  updator_type updator;
  return updator;
}

}  // namespace IndexedGridJsonHelper
}  // namespace Acts

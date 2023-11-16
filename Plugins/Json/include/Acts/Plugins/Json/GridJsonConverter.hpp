// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/detail/AxisFwd.hpp"

#include <iostream>

// Custom Json encoder/decoders.
namespace Acts {

namespace AxisJsonConverter {

/// Convert an axis to json
//
/// @param ia the axis
///
/// @return a json object to represent the axis
nlohmann::json toJson(const IAxis& ia);

}  // namespace AxisJsonConverter

namespace GridJsonConverter {

/// @brief Templated grid conversion to json
///
/// @tparam grid_type the type of the grid
/// @param grid the grid object
///
/// @return a json object to represent the grid
template <typename grid_type>
nlohmann::json toJson(const grid_type& grid) {
  nlohmann::json jGrid;

  auto axes = grid.axes();
  nlohmann::json jAxes;
  nlohmann::json jData;
  // 1D connections
  if constexpr (grid_type::DIM == 1u) {
    jAxes.push_back(AxisJsonConverter::toJson(*axes[0u]));
    for (unsigned int ib0 = 1u; ib0 <= axes[0u]->getNBins(); ++ib0) {
      typename grid_type::index_t lbin;
      lbin[0u] = ib0;
      jData.push_back(std::tie(lbin, grid.atLocalBins(lbin)));
    }
  }
  // 2D connections
  if constexpr (grid_type::DIM == 2u) {
    jAxes.push_back(AxisJsonConverter::toJson(*axes[0u]));
    jAxes.push_back(AxisJsonConverter::toJson(*axes[1u]));

    for (unsigned int ib0 = 1u; ib0 <= axes[0u]->getNBins(); ++ib0) {
      for (unsigned int ib1 = 1u; ib1 <= axes[1u]->getNBins(); ++ib1) {
        typename grid_type::index_t lbin;
        lbin[0u] = ib0;
        lbin[1u] = ib1;
        jData.push_back(std::tie(lbin, grid.atLocalBins(lbin)));
      }
    }
  }
  jGrid["axes"] = jAxes;
  jGrid["data"] = jData;

  return jGrid;
}

/// @brief Templated grid conversion from json
///
/// @tparam axis_generator_type the type of the axis generator
///         which determines the grid type
///
/// @param jGrid the json object to represent the grid
/// @param aGenerator the axis generator
///
/// @note the axis generator also defines the grid dimension
///
/// @return a grid object
template <typename axis_generator_type>
auto fromJson(const nlohmann::json& jGrid,
              const axis_generator_type& aGenerator) {
  // Generate the grid
  using GridType =
      typename axis_generator_type::template grid_type<std::vector<size_t>>;
  GridType grid(aGenerator());
  nlohmann::json jData = jGrid["data"];
  // Index filling
  if constexpr (GridType::DIM == 1u) {
    for (const auto& jd : jData) {
      std::array<size_t, 1u> lbin = jd[0u];
      std::vector<size_t> values = jd[1u];
      grid.atLocalBins(lbin) = values;
    }
  }
  if constexpr (GridType::DIM == 2u) {
    for (const auto& jd : jData) {
      std::array<size_t, 2u> lbin = jd[0u];
      std::vector<size_t> values = jd[1u];
      grid.atLocalBins(lbin) = values;
    }
  }

  return grid;
}

}  // namespace GridJsonConverter

/// @cond
NLOHMANN_JSON_SERIALIZE_ENUM(Acts::detail::AxisBoundaryType,
                             {{Acts::detail::AxisBoundaryType::Bound, "Bound"},
                              {Acts::detail::AxisBoundaryType::Open, "Open"},
                              {Acts::detail::AxisBoundaryType::Closed,
                               "Closed"}})

NLOHMANN_JSON_SERIALIZE_ENUM(Acts::detail::AxisType,
                             {{Acts::detail::AxisType::Equidistant,
                               "Equidistant"},
                              {Acts::detail::AxisType::Variable, "Variable"}})

/// @endcond

}  // namespace Acts

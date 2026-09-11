// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/AnyGridView.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/IGrid.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"
#include "ActsPlugins/Json/TrackParametersJsonConverter.hpp"

#include <iostream>
#include <stdexcept>

namespace Acts {

/// @addtogroup json_plugin
/// @{

/// @cond
NLOHMANN_JSON_SERIALIZE_ENUM(Acts::AxisBoundaryType,
                             {{Acts::AxisBoundaryType::Bound, "Bound"},
                              {Acts::AxisBoundaryType::Open, "Open"},
                              {Acts::AxisBoundaryType::Closed, "Closed"}})

NLOHMANN_JSON_SERIALIZE_ENUM(Acts::AxisType,
                             {{Acts::AxisType::Equidistant, "Equidistant"},
                              {Acts::AxisType::Variable, "Variable"}})
/// @endcond

namespace AxisJsonConverter {

/// Convert an axis to json
//
/// @param ia the axis
///
/// @return a json object to represent the axis
nlohmann::json toJson(const IAxis& ia);

/// Convert an axis from json
//
/// @param jAxis the serialized axis
///
/// @return axis object pointer
std::unique_ptr<Acts::IAxis> fromJson(const nlohmann::json& jAxis);

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
  for (unsigned int ia = 0u; ia < grid_type::DIM; ++ia) {
    auto jAxis = AxisJsonConverter::toJson(*axes[ia]);
    jAxes.push_back(jAxis);
  }
  jGrid["axes"] = jAxes;

  nlohmann::json jData;
  // 1D connections
  if constexpr (grid_type::DIM == 1u) {
    for (unsigned int ib0 = 1u; ib0 <= axes[0u]->getNBins(); ++ib0) {
      typename grid_type::index_t lbin;
      lbin[0u] = ib0;
      jData.push_back(std::tie(lbin, grid.atLocalBins(lbin)));
    }
  }
  // 2D connections
  if constexpr (grid_type::DIM == 2u) {
    for (unsigned int ib0 = 1u; ib0 <= axes[0u]->getNBins(); ++ib0) {
      for (unsigned int ib1 = 1u; ib1 <= axes[1u]->getNBins(); ++ib1) {
        typename grid_type::index_t lbin;
        lbin[0u] = ib0;
        lbin[1u] = ib1;
        jData.push_back(std::tie(lbin, grid.atLocalBins(lbin)));
      }
    }
  }
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
template <typename axis_generator_type,
          typename value_type = std::vector<std::size_t>>
auto fromJson(const nlohmann::json& jGrid,
              const axis_generator_type& aGenerator) {
  // Generate the grid
  using GridType = typename axis_generator_type::template grid_type<value_type>;
  GridType grid(aGenerator());
  nlohmann::json jData = jGrid["data"];
  // Index filling
  if constexpr (GridType::DIM == 1u) {
    for (const auto& jd : jData) {
      std::array<std::size_t, 1u> lbin = jd[0u];
      if (!jd[1u].is_null()) {
        grid.atLocalBins(lbin) = jd[1u].get<value_type>();
      }
    }
  }
  if constexpr (GridType::DIM == 2u) {
    for (const auto& jd : jData) {
      std::array<std::size_t, 2u> lbin = jd[0u];
      if (!jd[1u].is_null()) {
        grid.atLocalBins(lbin) = jd[1u].get<value_type>();
      }
    }
  }
  return grid;
}

/// @brief Type-erased grid conversion to json
///
/// Unlike @c toJson(const grid_type&), this works from the type-erased
/// @c IGrid / @c AnyGridConstView interfaces, so it does not need the
/// concrete axis types of the grid to be known at compile time.
///
/// @tparam value_type the type of the grid payload
/// @param grid the type-erased grid
/// @param view type-erased const view onto the grid values
///
/// @return a json object to represent the grid
template <typename value_type>
nlohmann::json toJsonAny(const IGrid& grid, AnyGridConstView<value_type> view) {
  std::size_t dim = grid.dimensions();
  if (dim != 1u && dim != 2u) {
    throw std::invalid_argument(
        "GridJsonConverter::toJsonAny: only 1D and 2D grids are supported.");
  }

  nlohmann::json jGrid;

  nlohmann::json jAxes;
  for (std::size_t ia = 0u; ia < dim; ++ia) {
    jAxes.push_back(AxisJsonConverter::toJson(grid.axis(ia)));
  }
  jGrid["axes"] = jAxes;

  nlohmann::json jData;
  if (dim == 1u) {
    for (std::size_t ib0 = 1u; ib0 <= grid.axis(0u).getNBins(); ++ib0) {
      IGrid::AnyIndexType lbin = {ib0};
      jData.push_back(std::tie(lbin, view.atLocalBins(lbin)));
    }
  } else {
    for (std::size_t ib0 = 1u; ib0 <= grid.axis(0u).getNBins(); ++ib0) {
      for (std::size_t ib1 = 1u; ib1 <= grid.axis(1u).getNBins(); ++ib1) {
        IGrid::AnyIndexType lbin = {ib0, ib1};
        jData.push_back(std::tie(lbin, view.atLocalBins(lbin)));
      }
    }
  }
  jGrid["data"] = jData;

  return jGrid;
}

}  // namespace GridJsonConverter

/// @}
}  // namespace Acts

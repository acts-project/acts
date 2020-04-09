// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialGridHelper.hpp, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include <stdexcept>

#include "Acts/Material/AccumulatedVolumeMaterial.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

namespace Acts {

using RecordedMaterialPoint =
    std::vector<std::pair<Acts::MaterialProperties, Acts::Vector3D>>;
using EAxis = Acts::detail::EquidistantAxis;
using Grid2D =
    Acts::detail::Grid<Acts::AccumulatedVolumeMaterial, EAxis, EAxis>;
using Grid3D =
    Acts::detail::Grid<Acts::AccumulatedVolumeMaterial, EAxis, EAxis, EAxis>;
using MaterialGrid2D = Acts::detail::Grid<Acts::ActsVectorF<5>, EAxis, EAxis>;
using MaterialGrid3D =
    Acts::detail::Grid<Acts::ActsVectorF<5>, EAxis, EAxis, EAxis>;

/// @brief Helper method that creates the cache grid for the mapping. This
/// grid allows the collection of material at a the anchor points.
///
/// @param [in] gridAxis1 Axis data
/// @param [in] gridAxis2 Axis data
/// @note The data of the axes is given in the std::array as {minimum value,
/// maximum value, number of bins}
///
/// @return The grid
Grid2D createGrid(std::array<double, 3> gridAxis1,
                  std::array<double, 3> gridAxis2);

/// @brief Helper method that creates the cache grid for the mapping. This
/// grid allows the collection of material at a the anchor points.
///
/// @param [in] gridAxis1 Axis data
/// @param [in] gridAxis2 Axis data
/// @param [in] gridAxis3 Axis data
/// @note The data of the axes is given in the std::array as {minimum value,
/// maximum value, number of bins}
///
/// @return The grid
Grid3D createGrid(std::array<double, 3> gridAxis1,
                  std::array<double, 3> gridAxis2,
                  std::array<double, 3> gridAxis3);

//
/// @brief Searcher for closest point in 2D coordinate
///
/// @param [in] matPos Position of the material
/// @param [in] grid Grid that is used for the look-up
///
/// @return Local grid point with the closest distance to @p matPos
Grid2D::index_t mapMaterial2D(const Acts::Vector3D& matPos, const Grid2D& grid);

/// @brief Searcher for closest point in 3D coordinate
///
/// @param [in] matPos Position of the material
/// @param [in] grid Grid that is used for the look-up
///
/// @return Local grid point with the closest distance to @p matPos
Grid3D::index_t mapMaterial3D(const Acts::Vector3D& matPos, const Grid3D& grid);

/// @brief Create a 3DGrid using a BinUtility.
/// Also determine the coresponding global to local transform and grid mapping
/// function
///
/// @param [in] BinUtility of the volume to be mapped
/// @param [in] Global to local transform to be updated.
///
/// @return the 3D grid
Grid3D createGrid(
    const BinUtility& bins,
    std::function<Acts::Vector3D(Acts::Vector3D)>& transfoGlobalToLocal);

/// @brief Concatenate a set of material at arbitrary space points on a set of
/// grid points and produces a grid containing the averaged material values.
///
/// @param [in] grid The material collecting grid
/// @param [in] mPoints The set of material at the space points
/// @param [in] transfoGlobalToLocal tranformation from local to local
/// coordinate
/// @param [in] matchToGridPoint Function that assigns the space point
/// of @p mPoints to the grid points by its local index
///
/// @return The average material grid decomposed into classification numbers
MaterialGrid2D mapMaterialPoints(
    Grid2D& grid, const Acts::RecordedMaterialPoint& mPoints,
    std::function<Acts::Vector3D(Acts::Vector3D)>& transfoGlobalToLocal,
    const std::function<Grid2D::index_t(const Acts::Vector3D&, const Grid2D&)>&
        matchToGridPoint);

/// @brief Concatenate a set of material at arbitrary space points on a set of
/// grid points and produces a grid containing the averaged material values.
///
/// @param [in] grid The material collecting grid
/// @param [in] mPoints The set of material at the space points
/// @param [in] transfoGlobalToLocal tranformation from local to local
/// coordinate
/// @param [in] matchToGridPoint Function that assigns the space point
/// of @p mPoints to the grid points by its local index
///
/// @return The average material grid decomposed into classification numbers
MaterialGrid3D mapMaterialPoints(
    Grid3D& grid, const Acts::RecordedMaterialPoint& mPoints,
    std::function<Acts::Vector3D(Acts::Vector3D)>& transfoGlobalToLocal,
    const std::function<Grid3D::index_t(const Acts::Vector3D&, const Grid3D&)>&
        matchToGridPoint);

}  // namespace Acts

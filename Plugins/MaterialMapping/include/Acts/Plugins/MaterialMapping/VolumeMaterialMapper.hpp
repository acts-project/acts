// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <vector>
#include "Acts/Material/Material.hpp"
#include "Acts/Plugins/MaterialMapping/AccumulatedVolumeMaterial.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

/// Convenience functions to ease creation of and Acts::InterpolatedMaterialMap
/// and to avoid code duplication. Currently implemented for the two most common
/// formats: rz and xyz.

namespace Acts {

/// @brief This function creates a discrete material map
///
/// @param [in] gridAxis1 Axis data
/// @param [in] gridAxis2 Axis data
/// @note The data of the axes is given in the std::array as {minimum value,
/// maximum value, number of bins}
/// @param [in] mPoints The set of material at the space points
/// @param [in] matchToGridPoint Function that assigns the space point
/// of @p mPoints to the grid points by its local index
///
/// @return The map
detail::Grid<ActsVectorF<5>, detail::EquidistantAxis, detail::EquidistantAxis>
createMaterialGrid(
    std::array<double, 3>                             gridAxis1,
    std::array<double, 3>                             gridAxis2,
    const std::vector<std::pair<Material, Vector3D>>& mPoints,
    const std::function<
        detail::Grid<ActsVectorF<5>,
                     detail::EquidistantAxis,
                     detail::EquidistantAxis>::
            index_t(const Vector3D&,
                    const detail::Grid<AccumulatedVolumeMaterial,
                                       detail::EquidistantAxis,
                                       detail::EquidistantAxis>&)>&
        matchToGridPoint);

/// @brief This function creates a discrete material map
///
/// @param [in] gridAxis1 Axis data
/// @param [in] gridAxis2 Axis data
/// @param [in] gridAxis3 Axis data
/// @note The data of the axes is given in the std::array as {minimum value,
/// maximum value, number of bins}
/// @param [in] mPoints The set of material at the space points
/// @param [in] matchToGridPoint Function that assigns the space point
/// of @p mPoints to the grid points by its local index
///
/// @return The map
detail::Grid<ActsVectorF<5>,
             detail::EquidistantAxis,
             detail::EquidistantAxis,
             detail::EquidistantAxis>
createMaterialGrid(
    std::array<double, 3>                             gridAxis1,
    std::array<double, 3>                             gridAxis2,
    std::array<double, 3>                             gridAxis3,
    const std::vector<std::pair<Material, Vector3D>>& mPoints,
    const std::function<
        detail::Grid<ActsVectorF<5>,
                     detail::EquidistantAxis,
                     detail::EquidistantAxis,
                     detail::EquidistantAxis>::
            index_t(const Vector3D&,
                    const detail::Grid<AccumulatedVolumeMaterial,
                                       detail::EquidistantAxis,
                                       detail::EquidistantAxis,
                                       detail::EquidistantAxis>&)>&
        matchToGridPoint);
}  // namespace Acts
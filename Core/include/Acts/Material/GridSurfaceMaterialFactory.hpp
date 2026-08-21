// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"

#include <memory>
#include <vector>

namespace Acts::GridSurfaceMaterialFactory {

/// Create a @c GridSurfaceMaterial with direct storage from two axes
///
/// @param axis0 the axis in direction 0
/// @param axis1 the axis in direction 1
/// @param payload the material payload, one slab per regular bin, column
///        major, i.e. [i0][i1]
/// @return a unique pointer to the created surface material
std::unique_ptr<GridSurfaceMaterial> createDirect(
    const IAxis& axis0, const IAxis& axis1,
    const std::vector<std::vector<MaterialSlab>>& payload);

/// Create a @c GridSurfaceMaterial with direct storage by resolving a
/// multi-axis spec against a surface
///
/// This follows the same pattern as @c ProtoGridSurfaceMaterial: the
/// (possibly deferred) axis specs are resolved against @p surface via
/// @c resolveMultiAxis, exactly as is done for a @c ProtoGridSurfaceMaterial
/// in @c BinnedSurfaceMaterialAccumulator. Binning restricted to a single
/// local direction is expressed by a single-bin spec in the other direction.
///
/// @param binning the 2D multi-axis binning spec (deferred axes are
///        resolved against @p surface)
/// @param surface the surface to resolve the deferred axes against
/// @param payload the material payload, one slab per regular bin, column
///        major, i.e. [i0][i1]
/// @return a unique pointer to the created surface material
std::unique_ptr<GridSurfaceMaterial> createDirect(
    const MultiAxisSpec2D& binning, const Surface& surface,
    const std::vector<std::vector<MaterialSlab>>& payload);

/// Create a @c GridSurfaceMaterial with locally indexed storage from two axes
///
/// @param axis0 the axis in direction 0
/// @param axis1 the axis in direction 1
/// @param material the locally owned material vector, addressed by @p payload
/// @param payload the index payload, one index per regular bin, column
///        major, i.e. [i0][i1]
/// @return a unique pointer to the created surface material
std::unique_ptr<GridSurfaceMaterial> createIndexed(
    const IAxis& axis0, const IAxis& axis1, std::vector<MaterialSlab> material,
    const std::vector<std::vector<std::size_t>>& payload);

/// Create a @c GridSurfaceMaterial with locally indexed storage from a
/// multi-axis spec resolved against a surface
///
/// @param binning the 2D multi-axis binning spec (deferred axes are
///        resolved against @p surface)
/// @param surface the surface to resolve the deferred axes against
/// @param material the locally owned material vector, addressed by @p payload
/// @param payload the index payload, one index per regular bin, column
///        major, i.e. [i0][i1]
/// @return a unique pointer to the created surface material
std::unique_ptr<GridSurfaceMaterial> createIndexed(
    const MultiAxisSpec2D& binning, const Surface& surface,
    std::vector<MaterialSlab> material,
    const std::vector<std::vector<std::size_t>>& payload);

/// Create a @c GridSurfaceMaterial with globally indexed storage from two axes
///
/// @param axis0 the axis in direction 0
/// @param axis1 the axis in direction 1
/// @param material the (possibly shared) globally owned material vector,
///        addressed by @p payload
/// @param payload the index payload, one index per regular bin, column
///        major, i.e. [i0][i1]
/// @return a unique pointer to the created surface material
std::unique_ptr<GridSurfaceMaterial> createGloballyIndexed(
    const IAxis& axis0, const IAxis& axis1,
    std::shared_ptr<std::vector<MaterialSlab>> material,
    const std::vector<std::vector<std::size_t>>& payload);

/// Create a @c GridSurfaceMaterial with globally indexed storage from a
/// multi-axis spec resolved against a surface
///
/// @param binning the binning specification for the grid
/// @param surface the surface to which the material is applied
/// @param material the (possibly shared) globally owned material vector,
///        addressed by @p payload
/// @param payload the index payload, one index per regular bin, column
///        major, i.e. [i0][i1]
/// @return a unique pointer to the created surface material
std::unique_ptr<GridSurfaceMaterial> createGloballyIndexed(
    const MultiAxisSpec2D& binning, const Surface& surface,
    std::shared_ptr<std::vector<MaterialSlab>> material,
    const std::vector<std::vector<std::size_t>>& payload);

}  // namespace Acts::GridSurfaceMaterialFactory

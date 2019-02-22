// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// VolumeMaterialMapper.hpp, Acts project MaterialPlugins
///////////////////////////////////////////////////////////////////

#pragma once

#include <vector>
#include "Acts/Material/Material.hpp"
#include "Acts/Plugins/MaterialMapping/AccumulatedVolumeMaterial.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

namespace Acts {

/// @brief This class serves to concatenate material evaluated at arbitrary
/// space points in a TrackingVolume on a provided set of grid points. The
/// procedure creates mean material values at these points and can be used to
/// create an interpolated grid of material.
/// @note Since this class just associates an arbitrary cloud of space points
/// with the material evaluated at these points to a user defined grid, this
/// procedure is a pure pre-processing of an applicable material getter. This
/// last step requires additional interpolation between the grid points.
class VolumeMaterialMapper
{
public:
  using RecordedMaterial = std::vector<std::pair<Material, Vector3D>>;
  using EAxis            = detail::EquidistantAxis;
  using Grid2D = detail::Grid<AccumulatedVolumeMaterial, EAxis, EAxis>;
  using Grid3D = detail::Grid<AccumulatedVolumeMaterial, EAxis, EAxis, EAxis>;
  using MaterialGrid2D = detail::Grid<ActsVectorF<5>, EAxis, EAxis>;
  using MaterialGrid3D = detail::Grid<ActsVectorF<5>, EAxis, EAxis, EAxis>;

  /// @brief Default constructor
  VolumeMaterialMapper() = default;

  /// @brief Default destructor
  ~VolumeMaterialMapper() = default;

  /// @brief Helper method that creates the cache grid for the mapping. This
  /// grid allows the collection of material at a the anchro points.
  ///
  /// @param [in] gridAxis1 Vector of grid points in the arbitrary first
  /// dimension
  /// @param [in] gridAxis2 Vector of grid points in the arbitrary second
  /// dimension
  ///
  /// @return The grid
  Grid2D
  createGrid(std::vector<double> gridAxis1,
             std::vector<double> gridAxis2) const;

  /// @brief Helper method that creates the cache grid for the mapping. This
  /// grid allows the collection of material at a the anchro points.
  ///
  /// @param [in] gridAxis1 Vector of grid points in the arbitrary first
  /// dimension
  /// @param [in] gridAxis2 Vector of grid points in the arbitrary second
  /// dimension
  /// @param [in] gridAxis3 Vector of grid points in the arbitrary third
  /// dimension
  ///
  /// @return The grid
  Grid3D
  createGrid(std::vector<double> gridAxis1,
             std::vector<double> gridAxis2,
             std::vector<double> gridAxis3) const;

  /// @brief Concatenate a set of material at arbitrary space points on a set of
  /// grid points and produces a grid containing the averaged material values.
  ///
  /// @param [in] grid The material collecting grid
  /// @param [in] mPoints The set of material at the space points
  /// @param [in] concatenateToGridPoints Function that assigns the space point
  /// of @p mPoints to the grid points by its local index
  ///
  /// @return The average material grid decomposed into classification numbers
  MaterialGrid2D
  mapMaterialPoints(
      Grid2D&                 grid,
      const RecordedMaterial& mPoints,
      const std::function<Grid2D::index_t(const Vector3D&, const Grid2D&)>&
          concatenateToGridPoints) const;

  /// @brief Concatenate a set of material at arbitrary space points on a set of
  /// grid points and produces a grid containing the averaged material values.
  ///
  /// @param [in] grid The material collecting grid
  /// @param [in] mPoints The set of material at the space points
  /// @param [in] concatenateToGridPoints Function that assigns the space point
  /// of @p mPoints to the grid points by its local index
  ///
  /// @return The average material grid decomposed into classification numbers
  MaterialGrid3D
  mapMaterialPoints(
      Grid3D&                 grid,
      const RecordedMaterial& mPoints,
      const std::function<Grid3D::index_t(const Vector3D&, const Grid3D&)>&
          concatenateToGridPoints) const;
};
}

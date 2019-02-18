// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <functional>
#include "Acts/Material/concept/AnyFieldLookup.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Interpolation.hpp"
#include "Acts/Utilities/concept/AnyGrid.hpp"
#include "Acts/Material/Material.hpp"

namespace Acts {


class InterpolatedMaterialMap final
{

  template <typename T, class Point, size_t DIM>
  using AnyGrid_t = concept::AnyNDimInterpGrid<T, Point, DIM>;

  template <unsigned int DIM_POS>
  struct FieldCell
  {
    static constexpr unsigned int N = 1 << DIM_POS;

  public:
    FieldCell(std::function<ActsVectorD<DIM_POS>(const Vector3D&)> transformPos,
              std::array<double, DIM_POS> lowerLeft,
              std::array<double, DIM_POS> upperRight,
              std::array<Vector3D, N>     fieldValues)
      : m_transformPos(std::move(transformPos))
      , m_lowerLeft(std::move(lowerLeft))
      , m_upperRight(std::move(upperRight))
      , m_fieldValues(std::move(fieldValues))
    {
    }

    Vector3D
    getMaterial(const Vector3D& position) const
    {
      // defined in Interpolation.hpp
      return interpolate(
          m_transformPos(position), m_lowerLeft, m_upperRight, m_fieldValues);
    }

    bool
    isInside(const Vector3D& position) const
    {
      const auto& gridCoordinates = m_transformPos(position);
      for (unsigned int i = 0; i < DIM_POS; ++i) {
        if (gridCoordinates[i] < m_lowerLeft.at(i)
            || gridCoordinates[i] >= m_upperRight.at(i)) {
          return false;
        }
      }
      return true;
    }

  private:
    /// geometric transformation applied to global 3D positions
    std::function<ActsVectorD<DIM_POS>(const Vector3D&)> m_transformPos;

    /// generalized lower-left corner of the confining hyper-box
    std::array<double, DIM_POS> m_lowerLeft;

    /// generalized upper-right corner of the confining hyper-box
    std::array<double, DIM_POS> m_upperRight;

    /// @brief magnetic field vectors at the hyper-box corners
    ///
    /// @note These values must be order according to the prescription detailed
    ///       in Acts::interpolate.
    std::array<Vector3D, N> m_fieldValues;
  };

public:
// DIM_POS: Local frame for parameter look-up
// DIM_BFIELD: Returning dimension of B-field after look-up
  template <unsigned int DIM_POS>
  struct FieldMapper
  {
  public:
    using Grid_t
        = AnyGrid_t<Material, ActsVectorD<DIM_POS>, DIM_POS>;
    
    FieldMapper(
        std::function<ActsVectorD<DIM_POS>(const Vector3D&)> transformPos,
        Grid_t grid)
      : m_transformPos(std::move(transformPos))
      , m_grid(std::move(grid))
    {
    }

    Material
    getMaterial(const Vector3D& position) const
    {
      return m_grid.interpolate(m_transformPos(position)); // TODO: interpolate must return the material - should be covered by the template parameter of the grid
    }

    FieldCell<DIM_POS>
    getMaterialCell(const Vector3D& position) const
    {
      const auto& gridPosition = m_transformPos(position);
      size_t      bin          = m_grid.getGlobalBinIndex(gridPosition);
      const auto& indices      = m_grid.getLocalBinIndices(bin);
      const auto& lowerLeft    = m_grid.getLowerLeftBinEdge(indices);
      const auto& upperRight   = m_grid.getUpperRightBinEdge(indices);

      // loop through all corner points
      constexpr size_t nCorners = 1 << DIM_POS;
      std::array<Vector3D, nCorners> neighbors;
      const auto& cornerIndices = m_grid.closestPointsIndices(gridPosition);

      size_t i = 0;
      for (size_t index : cornerIndices) {
        neighbors.at(i++) = m_grid.at(index); // TODO: this might work
      }

      return FieldCell<DIM_POS>(
          m_transformPos, lowerLeft, upperRight, std::move(neighbors));
    }

    std::vector<size_t>
    getNBins() const
    {
      auto nBinsArray = m_grid.getNBins();
      return std::vector<size_t>(nBinsArray.begin(), nBinsArray.end());
    }
// TODO: min/max should be working, but this might be checked
    std::vector<double>
    getMin() const
    {
      auto minArray = m_grid.getMin();
      return std::vector<double>(minArray.begin(), minArray.end());
    }

    std::vector<double>
    getMax() const
    {
      auto maxArray = m_grid.getMax();
      return std::vector<double>(maxArray.begin(), maxArray.end());
    }

    bool
    isInside(const Vector3D& position) const
    {
      return m_grid.isInside(m_transformPos(position));
    }

    const Grid_t&
    getGrid() const
    {
      return m_grid;
    }

  private:
    /// geometric transformation applied to global 3D positions
    std::function<ActsVectorD<DIM_POS>(const Vector3D&)> m_transformPos;
    /// grid storing magnetic field values
    Grid_t m_grid;
  };

  /// @brief configuration object for magnetic field interpolation
  struct Config
  {
// TODO: might be useless
    //~ /// @brief global B-field scaling factor
    //~ ///
    //~ /// @note Negative values for @p scale are accepted and will invert the
    //~ ///       direction of the magnetic field.
    //~ double scale = 1.;

    /// @brief object for global coordinate transformation and interpolation
    ///
    /// This object performs the mapping of the global 3D coordinates onto the
    /// field grid and the interpolation of the field values on close-by grid
    /// points.
    concept::AnyFieldLookup<> mapper;
  };

  struct Cache
  {
    concept::AnyFieldCell<> fieldCell;
    bool                    initialized = false;
  };

  InterpolatedMaterialMap(Config config) : m_config(std::move(config)) {}

  Config
  getConfiguration() const
  {
    return m_config;
  }

  Vector3D
  getMaterial(const Vector3D& position) const
  {
    return m_config.mapper.getMaterial(position);
  }

  Vector3D
  getMaterial(const Vector3D& position, Cache& cache) const
  {
    if (!cache.initialized || !cache.fieldCell.isInside(position)) {
      cache.fieldCell   = getMaterialCell(position);
      cache.initialized = true;
    }
    return cache.fieldCell.getMaterial(position);
  }

  Vector3D
  getMaterialGradient(const Vector3D& position,
                   ActsMatrixD<3, 3>& /*derivative*/) const
  {
    return m_config.mapper.getMaterial(position);
  }

  Vector3D
  getMaterialGradient(const Vector3D& position,
                   ActsMatrixD<3, 3>& /*derivative*/,
                   Cache& /*cache*/) const
  {
    return m_config.mapper.getMaterial(position);
  }

  concept::AnyFieldLookup<>
  getMapper() const
  {
    return m_config.mapper;
  }

  bool
  isInside(const Vector3D& position) const
  {
    return m_config.mapper.isInside(position);
  }

  void
  setConfiguration(const Config& config)
  {
    m_config = config;
  }

private:
  concept::AnyFieldCell<>
  getMaterialCell(const Vector3D& position) const
  {
    return m_config.mapper.getMaterialCell(position);
  }

  Config m_config;
};

}  // namespace Acts

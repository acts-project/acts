// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/BFieldProvider.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Utilities/Interpolation.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

#include <functional>
#include <optional>
#include <vector>

namespace Acts {

/// @brief struct for mapping global 3D positions to field values
///
/// @tparam DIM_POS Dimensionality of position in magnetic field map
/// @tparam DIM_BFIELD Dimensionality of BField in magnetic field map
///
/// Global 3D positions are transformed into a @c DIM_POS Dimensional
/// vector which
/// is used to look up the magnetic field value in the underlying field map.
template <typename G>
struct InterpolatedBFieldMapper {
 public:
  using Grid_t = G;
  using FieldType = typename Grid_t::value_type;
  static constexpr size_t DIM_POS = Grid_t::DIM;

  /// @brief struct representing smallest grid unit in magnetic field grid
  ///
  /// This type encapsulate all required information to perform linear
  /// interpolation of magnetic field values within a confined 3D volume.
  struct FieldCell {
    /// number of corner points defining the confining hyper-box
    static constexpr unsigned int N = 1 << DIM_POS;

   public:
    /// @brief default constructor
    ///
    /// @param [in] transform   mapping of global 3D coordinates onto grid space
    /// @param [in] lowerLeft   generalized lower-left corner of hyper box
    ///                         (containing the minima of the hyper box along
    ///                         each Dimension)
    /// @param [in] upperRight  generalized upper-right corner of hyper box
    ///                         (containing the maxima of the hyper box along
    ///                         each Dimension)
    /// @param [in] fieldValues field values at the hyper box corners sorted in
    ///                         the canonical order defined in Acts::interpolate
    FieldCell(std::function<ActsVector<DIM_POS>(const Vector3&)> transformPos,
              std::array<double, DIM_POS> lowerLeft,
              std::array<double, DIM_POS> upperRight,
              std::array<Vector3, N> fieldValues)
        : m_transformPos(std::move(transformPos)),
          m_lowerLeft(std::move(lowerLeft)),
          m_upperRight(std::move(upperRight)),
          m_fieldValues(std::move(fieldValues)) {}

    /// @brief retrieve field at given position
    ///
    /// @param [in] position global 3D position
    /// @return magnetic field value at the given position
    ///
    /// @pre The given @c position must lie within the current field cell.
    Vector3 getField(const Vector3& position) const {
      // defined in Interpolation.hpp
      return interpolate(m_transformPos(position), m_lowerLeft, m_upperRight,
                         m_fieldValues);
    }

    /// @brief check whether given 3D position is inside this field cell
    ///
    /// @param [in] position global 3D position
    /// @return @c true if position is inside the current field cell,
    ///         otherwise @c false
    bool isInside(const Vector3& position) const {
      const auto& gridCoordinates = m_transformPos(position);
      for (unsigned int i = 0; i < DIM_POS; ++i) {
        if (gridCoordinates[i] < m_lowerLeft.at(i) ||
            gridCoordinates[i] >= m_upperRight.at(i)) {
          return false;
        }
      }
      return true;
    }

   private:
    /// geometric transformation applied to global 3D positions
    std::function<ActsVector<DIM_POS>(const Vector3&)> m_transformPos;

    /// generalized lower-left corner of the confining hyper-box
    std::array<double, DIM_POS> m_lowerLeft;

    /// generalized upper-right corner of the confining hyper-box
    std::array<double, DIM_POS> m_upperRight;

    /// @brief magnetic field vectors at the hyper-box corners
    ///
    /// @note These values must be order according to the prescription detailed
    ///       in Acts::interpolate.
    std::array<Vector3, N> m_fieldValues;
  };

  /// @brief default constructor
  ///
  /// @param [in] transformPos mapping of global 3D coordinates (cartesian)
  /// onto grid space
  /// @param [in] transformBField calculating the global 3D coordinates
  /// (cartesian) of the magnetic field with the local n dimensional field and
  /// the global 3D position as input
  /// @param [in] grid      grid storing magnetic field values
  InterpolatedBFieldMapper(
      std::function<ActsVector<DIM_POS>(const Vector3&)> transformPos,
      std::function<Vector3(const FieldType&, const Vector3&)> transformBField,
      Grid_t grid)
      : m_transformPos(std::move(transformPos)),
        m_transformBField(std::move(transformBField)),
        m_grid(std::move(grid)) {}

  /// @brief retrieve field at given position
  ///
  /// @param [in] position global 3D position
  /// @return magnetic field value at the given position
  ///
  /// @pre The given @c position must lie within the range of the underlying
  ///      magnetic field map.
  Vector3 getField(const Vector3& position) const {
    return m_transformBField(m_grid.interpolate(m_transformPos(position)),
                             position);
  }

  /// @brief retrieve field cell for given position
  ///
  /// @param [in] position global 3D position
  /// @return field cell containing the given global position
  ///
  /// @pre The given @c position must lie within the range of the underlying
  ///      magnetic field map.
  FieldCell getFieldCell(const Vector3& position) const {
    const auto& gridPosition = m_transformPos(position);
    const auto& indices = m_grid.localBinsFromPosition(gridPosition);
    const auto& lowerLeft = m_grid.lowerLeftBinEdge(indices);
    const auto& upperRight = m_grid.upperRightBinEdge(indices);

    // loop through all corner points
    constexpr size_t nCorners = 1 << DIM_POS;
    std::array<Vector3, nCorners> neighbors;
    const auto& cornerIndices = m_grid.closestPointsIndices(gridPosition);

    size_t i = 0;
    for (size_t index : cornerIndices) {
      neighbors.at(i++) = m_transformBField(m_grid.at(index), position);
    }

    return FieldCell(m_transformPos, lowerLeft, upperRight,
                     std::move(neighbors));
  }

  /// @brief get the number of bins for all axes of the field map
  ///
  /// @return vector returning number of bins for all field map axes
  std::vector<size_t> getNBins() const {
    auto nBinsArray = m_grid.numLocalBins();
    return std::vector<size_t>(nBinsArray.begin(), nBinsArray.end());
  }

  /// @brief get the minimum value of all axes of the field map
  ///
  /// @return vector returning the minima of all field map axes
  std::vector<double> getMin() const {
    auto minArray = m_grid.minPosition();
    return std::vector<double>(minArray.begin(), minArray.end());
  }

  /// @brief get the maximum value of all axes of the field map
  ///
  /// @return vector returning the maxima of all field map axes
  std::vector<double> getMax() const {
    auto maxArray = m_grid.maxPosition();
    return std::vector<double>(maxArray.begin(), maxArray.end());
  }

  /// @brief check whether given 3D position is inside look-up domain
  ///
  /// @param [in] position global 3D position
  /// @return @c true if position is inside the defined look-up grid,
  ///         otherwise @c false
  bool isInside(const Vector3& position) const {
    return m_grid.isInside(m_transformPos(position));
  }

  /// @brief Get a const reference on the underlying grid structure
  ///
  /// @return grid reference
  const Grid_t& getGrid() const { return m_grid; }

 private:
  /// geometric transformation applied to global 3D positions
  std::function<ActsVector<DIM_POS>(const Vector3&)> m_transformPos;
  /// Transformation calculating the global 3D coordinates (cartesian) of the
  /// magnetic field with the local n dimensional field and the global 3D
  /// position as input
  std::function<Vector3(const FieldType&, const Vector3&)> m_transformBField;
  /// grid storing magnetic field values
  Grid_t m_grid;
};

/// @ingroup MagneticField
/// @brief interpolate magnetic field value from field values on a given grid
///
/// This class implements a magnetic field service which is initialized by a
/// field map defined by:
/// - a list of field values on a regular grid in some n-Dimensional space,
/// - a transformation of global 3D coordinates onto this n-Dimensional
/// space.
/// - a transformation of local n-Dimensional magnetic field coordinates into
/// global (cartesian) 3D coordinates
///
/// The magnetic field value for a given global position is then determined by:
/// - mapping the position onto the grid,
/// - looking up the magnetic field values on the closest grid points,
/// - doing a linear interpolation of these magnetic field values.
template <typename Mapper_t>
class InterpolatedBFieldMap final : public BFieldProvider {
 public:
  /// @brief configuration object for magnetic field interpolation
  struct Config {
    Config(Mapper_t m) : mapper(m) {}

    /// @brief global B-field scaling factor
    ///
    /// @note Negative values for @p scale are accepted and will invert the
    ///       direction of the magnetic field.
    double scale = 1.;

    /// @brief object for global coordinate transformation and interpolation
    ///
    /// This object performs the mapping of the global 3D coordinates onto the
    /// field grid and the interpolation of the field values on close-by grid
    /// points.
    Mapper_t mapper;
  };

  struct Cache {
    /// @brief Constructor with magnetic field context
    ///
    /// @param mcfg the magnetic field context
    Cache(const MagneticFieldContext& /*mcfg*/) {}

    std::optional<typename Mapper_t::FieldCell> fieldCell;
    bool initialized = false;
  };

  /// @brief create interpolated magnetic field map
  ///
  /// @param [in] config configuration object
  InterpolatedBFieldMap(Config config) : m_config(std::move(config)) {}

  /// @brief get configuration object
  ///
  /// @return copy of the internal configuration object
  Config getConfiguration() const { return m_config; }

  /// @brief get global scaling factor for magnetic field
  ///
  /// @return global factor for scaling the magnetic field
  double getScale() const { return m_config.scale; }

  /// @brief convenience method to access underlying field mapper
  ///
  /// @return the field mapper
  Mapper_t getMapper() const { return m_config.mapper; }

  /// @brief check whether given 3D position is inside look-up domain
  ///
  /// @param [in] position global 3D position
  /// @return @c true if position is inside the defined BField map,
  ///         otherwise @c false
  bool isInside(const Vector3& position) const {
    return m_config.mapper.isInside(position);
  }

  /// @brief update configuration
  ///
  /// @param [in] config new configuration object
  void setConfiguration(const Config& config) { m_config = config; }

 public:
  /// @copydoc BFieldBase::makeCache(const MagneticFieldContext&)
  BFieldProvider::Cache makeCache(
      const MagneticFieldContext& mctx) const override {
    return BFieldProvider::Cache::make<Cache>(mctx);
  }

  /// @copydoc BFieldBase::getField(const Vector3&)
  Vector3 getField(const Vector3& position) const override {
    return m_config.mapper.getField(position);
  }

  /// @copydoc BFieldBase::getField(const Vector3&,BFieldBase::Cache&)
  Vector3 getField(const Vector3& position,
                   BFieldProvider::Cache& gcache) const override {
    Cache& cache = gcache.get<Cache>();
    if (!cache.fieldCell || !(*cache.fieldCell).isInside(position)) {
      cache.fieldCell = getFieldCell(position);
    }
    return (*cache.fieldCell).getField(position);
  }

  /// @copydoc BFieldBase::getFieldGradient(const Vector3&,ActsMatrix<3,3>&)
  ///
  /// @note currently the derivative is not calculated
  /// @todo return derivative
  Vector3 getFieldGradient(const Vector3& position,
                           ActsMatrix<3, 3>& /*derivative*/) const override {
    return m_config.mapper.getField(position);
  }

  /// @copydoc BFieldBase::getFieldGradient(const
  /// Vector3&,ActsMatrix<3,3>&,BFieldBase::Cache&)
  ///
  /// @note currently the derivative is not calculated
  /// @note Cache is not used currently
  /// @todo return derivative
  Vector3 getFieldGradient(const Vector3& position,
                           ActsMatrix<3, 3>& /*derivative*/,
                           BFieldProvider::Cache& /*cache*/) const override {
    return m_config.mapper.getField(position);
  }

 private:
  /// @brief retrieve field cell for given position
  ///
  /// @param [in] position global 3D position
  /// @return field cell containing the given global position
  ///
  /// @pre The given @c position must lie within the range of the underlying
  ///      magnetic field map.
  typename Mapper_t::FieldCell getFieldCell(const Vector3& position) const {
    return m_config.mapper.getFieldCell(position);
  }

  /// @brief configuration object
  Config m_config;
};

}  // namespace Acts

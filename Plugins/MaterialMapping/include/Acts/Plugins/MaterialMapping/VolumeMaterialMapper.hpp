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

#include "Acts/Material/Material.hpp"
#include "Acts/Plugins/MaterialMapping/AccumulatedVolumeMaterial.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

class TrackingVolume;

/// @brief This class serves to concatenate material evaluated at arbitrary
/// space points in a TrackingVolume on a provided set of grid points. The
/// procedure creates mean material values at these points and can be used to
/// create an interpolated grid of material.
class VolumeMaterialMapper
{
public:
  using RecordedMaterial = std::vector<std::pair<Material, Vector3D>>

      /// @struct State
      ///
      /// State struct which is used for the mapping prococess
      struct State
  {
    /// Storage of the grid points in each dimension. The structure is
    /// gridPointsPerAxis[dimension][index]
    std::vector<std::vector<double>> gridPointsPerAxis;
    /// Storage of the accumulated material obtained at each of the grid points.
    /// The vector has the length of all combinatorial grid points.
    std::vector<AccumulatedVolumeMaterial> accumulatedMaterial;
  };

  /// Delete the Default constructor
  VolumeMaterialMapper() = delete;

  /// Constructor with config object
  ///
  /// @param [in] log The logger
  VolumeMaterialMapper(std::unique_ptr<const Logger> slogger
                       = getDefaultLogger("SurfaceMaterialMapper",
                                          Logging::INFO));

  /// @brief Helper method that creates the cache for the mapping.
  ///
  /// @param [in] gridAxis1 Vector of grid points in the arbitrary first
  /// dimension
  /// @param [in] gridAxis2 Vector of grid points in the arbitrary second
  /// dimension
  /// @param [in] gridAxis3 Vector of grid points in the arbitrary third
  /// dimension
  ///
  /// @note The number of filled vectors describe the dimension inside the space
  /// of the latter look-up grid of the material. The dimension of the mapping
  /// is not fix and can be 1-3 dimensional depending of the amount of non-empty
  /// vectors.
  /// @note At this point there does not exist a connection of the entries of
  /// the vector to an axis. This connection becomes important in the
  /// VolumeMaterialMapper::mapMaterialPoints() function. Also the values in the
  /// vectors are sorted inside the function.
  /// @return State object
  State
  createState(const vector<double>& gridAxis1,
              const vector<double>& gridAxis2 = {},
              const vector<double>& gridAxis3 = {}) const;

  /// @brief Concatenate a set of material at arbitrary space points on a set of
  /// grid points
  ///
  /// @param [in] mState The current state map
  /// @param [in] mPoints The set of material at the space points
  /// @param [in] concatenateToGridPoints Function that assigns the space point
  /// of @p mPoints to the grid points
  void
  mapMaterialPoints(
      State&                  mState,
      const RecordedMaterial& mPoints,
      const std::function<unsigned int(const Vector3D&, const State&)>&
          concatenateToGridPoints) const;

  /// @brief Method to finalize the maps
  ///
  /// It calls the final run averaging and then transforms
  /// the AccumulatedSurface material class to a vector of materials with the
  /// same ordering as the mState.accumulatedMaterial vector.
  ///
  /// @param [in] mState State object
  ///
  /// @return Vector of material at each grid point
  std::vector<Material>
  finalizeMaps(State& mState) const;

private:
  /// The logging instance
  std::unique_ptr<const Logger> m_logger;
};
}

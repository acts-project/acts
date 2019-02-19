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

#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Plugins/MaterialMapping/AccumulatedVolumeMaterial.hpp"
#include "Acts/Plugins/MaterialMapping/RecordedMaterialTrack.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Material/Material.hpp"

namespace Acts {

class TrackingVolume;

class VolumeMaterialMapper
{
public:
  using StraightLinePropagator = Propagator<StraightLineStepper, Navigator>;

  /// @struct State
  ///
  /// Nested State struct which is used for the mapping prococess
  struct State
  {
	  // TODO: Require a functional description for association of space point with anchor point
	  std::vector<std::vector<double>> edgesPerAxis;
	  
	  std::vector<AccumulatedVolumeMaterial> accumulatedMaterial;
	  
    //~ /// The accumulated material per geometry ID
    //~ std::map<GeometryID, AccumulatedVolumeMaterial> accumulatedMaterial; // TODO: GeoID should be something related to an axis construct
    //~ /// The created surface material from it
    //~ std::map<GeometryID, std::unique_ptr<const SurfaceMaterial>> // TODO: must become a look-up of the grid points or since this keeps the binning be removed
        //~ surfaceMaterial;
  };

  /// Delete the Default constructor
  VolumeMaterialMapper() = delete;

  /// Constructor with config object
  ///
  /// @param cfg Configuration struct
  /// @param propagator The straight line propagator
  /// @param log The logger
  VolumeMaterialMapper(bool mapperDebugOutput                 dpg,
                        std::unique_ptr<const Logger> slogger
                        = getDefaultLogger("SurfaceMaterialMapper",
                                           Logging::INFO));

  /// @brief helper method that creates the cache for the mapping
  ///
  /// @param[in] tGeometry The geometry which should be mapped
  ///
  /// This method takes a TrackingGeometry,
  /// finds all surfaces with material proxis
  /// and returns you a Cache object tO be used
  // TODO: Sorting inside of function
  State
  createState(const vector<double>& edgeAxis1, const vector<double>& edgeAxis2 = {}, const vector<double>& edgeAxis3 = {}) const;

  /// @brief Method to finalize the maps
  ///
  /// It calls the final run averaging and then transforms
  /// the AccumulatedSurface material class to a surface material
  /// class type
  ///
  /// @param mState
  std::vector<Material>
  finalizeMaps(State& mState) const;

  /// Process/map a single track
  ///
  /// @param mState The current state map
  /// @param mTrack The material track to be mapped
  ///
  /// @note the RecordedMaterialProperties of the track are assumed
  /// to be ordered from the starting position along the starting direction
  // TODO: mTrack should become vector<material, vector3d>
  void
  mapMaterialTrack(State& mState, const RecordedMaterialTrack& mTrack, const std::function<unsigned int(const Vector3D&, const State&)>& concatenateToEdge) const;

private:

    /// Mapping output to debug stream
    bool m_mapperDebugOutput = false;

  /// The logging instance
  std::unique_ptr<const Logger> m_logger;
};
}

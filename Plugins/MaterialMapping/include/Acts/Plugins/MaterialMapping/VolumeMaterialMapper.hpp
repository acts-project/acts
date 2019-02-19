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
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/detail/Axis.hpp"

namespace Acts {

class TrackingVolume;

class VolumeMaterialMapper
{
public:
  using StraightLinePropagator = Propagator<StraightLineStepper, Navigator>;

  struct AssignedMaterialProperties
  {
    /// GeometryID of object were recorded material properties is assigned
    GeometryID geoID{0}; // TODO: should be replaced or removed (latter if the mapping is performed per volume)
    /// Projected position of the assigned material
    Vector3D assignedPosition{0., 0., 0.};
    /// The material information that is assigned (can be multiple steps)
    std::vector<RecordedMaterialProperties> assignedProperties; // TODO: The position stored in the values should just be mapped
    /// The incident angle correction due to particle incident
    double pathCorrection{0.}; // TODO: not sure if this can remain

    /// @brief Constructor for AssignedMaterialProperties
    ///
    /// @param gid The GeometryID of the surface
    /// @param pos The position of the assignment/intersection
    /// @param pc The pathCorrection to be applied (inverse)
    AssignedMaterialProperties(GeometryID gid, Vector3D pos, double pc)
      : geoID(gid), assignedPosition(std::move(pos)), pathCorrection(pc)
    {
    }
  };

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
  State
  createState(const vector<double>& edgeAxis1, const vector<double>& edgeAxis2 = {}, const vector<double>& edgeAxis3 = {}) const;

  /// @brief Method to finalize the maps
  ///
  /// It calls the final run averaging and then transforms
  /// the AccumulatedSurface material class to a surface material
  /// class type
  ///
  /// @param mState
  void
  finalizeMaps(State& mState) const;

  /// Process/map a single track
  ///
  /// @param mState The current state map
  /// @param mTrack The material track to be mapped
  ///
  /// @note the RecordedMaterialProperties of the track are assumed
  /// to be ordered from the starting position along the starting direction
  void
  mapMaterialTrack(State& mState, const RecordedMaterialTrack& mTrack) const;

private:

    /// Mapping output to debug stream
    bool m_mapperDebugOutput = false;

  /// The straight line propagator
  StraightLinePropagator m_propagator;

  /// The logging instance
  std::unique_ptr<const Logger> m_logger;
};
}

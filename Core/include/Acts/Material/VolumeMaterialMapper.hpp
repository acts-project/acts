// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/AccumulatedVolumeMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/VolumeCollector.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"
#include "Acts/Utilities/Logger.hpp"

/// Convenience functions to ease creation of and Acts::InterpolatedMaterialMap
/// and to avoid code duplication. Currently implemented for the two most common
/// formats: rz and xyz.

namespace Acts {

using RecordedMaterialPoint = std::vector<std::pair<Acts::MaterialProperties, Acts::Vector3D>>;

class VolumeMaterialMapper {
 public:
   using StraightLinePropagator = Propagator<StraightLineStepper, Navigator>;

   /// @brief selector for finding
   struct MaterialVolume {
     bool operator()(const TrackingVolume& vf) const {
       return (vf.volumeMaterial() != nullptr);
     }
   };

  /// @struct Config
  ///
  /// Nested Configuration struct for the material mapper
  struct Config {
    /// Mapping range
    std::array<double, 2> etaRange = {{-6., 6.}};
    /// Correct for empty bins (recommended)
    bool emptyBinCorrection = true;
    /// Mapping output to debug stream
    bool mapperDebugOutput = false;
  };

  /// @struct State
  ///
  /// Nested State struct which is used for the mapping prococess
  struct State {
    /// Constructor of the Sate with contexts
    State(std::reference_wrapper<const GeometryContext> gctx,
          std::reference_wrapper<const MagneticFieldContext> mctx)
        : geoContext(gctx), magFieldContext(mctx) {}

    /// The accumulated material per geometry ID
    std::map<GeometryID, RecordedMaterialPoint> recordedMaterial;

    /// The binning per geometry ID
    std::map<GeometryID, BinUtility> materialBin;

    /// The created surface material from it
    std::map<GeometryID, std::unique_ptr<const IVolumeMaterial>>
        volumeMaterial;

    /// The created surface material from it
    std::map<GeometryID, std::shared_ptr<const ISurfaceMaterial>>
        surfaceMaterial;

    /// Reference to the geometry context for the mapping
    std::reference_wrapper<const GeometryContext> geoContext;

    /// Reference to the magnetic field context
    std::reference_wrapper<const MagneticFieldContext> magFieldContext;

  };

  /// Delete the Default constructor
  VolumeMaterialMapper() = delete;

  /// Constructor with config object
  ///
  /// @param cfg Configuration struct
  /// @param propagator The straight line propagator
  /// @param log The logger
  VolumeMaterialMapper(const Config& cfg, StraightLinePropagator propagator,
                       std::unique_ptr<const Logger> slogger =
                            getDefaultLogger("VolumeMaterialMapper",
                                             Logging::INFO));

  /// @brief helper method that creates the cache for the mapping
  ///
  /// @param[in] tGeometry The geometry which should be mapped
  ///
  /// This method takes a TrackingGeometry,
  /// finds all surfaces with material proxis
  /// and returns you a Cache object tO be used
  State createState(const GeometryContext& gctx,
                    const MagneticFieldContext& mctx,
                    const TrackingGeometry& tGeometry) const;

  /// @brief Method to finalize the maps
  ///
  /// It calls the final run averaging and then transforms
  /// the AccumulatedSurface material class to a surface material
  /// class type
  ///
  /// @param mState
  void finalizeMaps(State& mState) const;

  /// Process/map a single track
  ///
  /// @param mState The current state map
  /// @param mTrack The material track to be mapped
  ///
  /// @note the RecordedMaterialProperties of the track are assumed
  /// to be ordered from the starting position along the starting direction
  void mapMaterialTrack(State& mState, RecordedMaterialTrack& mTrack) const;

 private:
  /// @brief finds all surfaces with ProtoSurfaceMaterial of a volume
  ///
  /// @param mState The state to be filled
  /// @param tVolume is current TrackingVolume
  void resolveMaterialVolume(State& mState,
                               const TrackingVolume& tVolume) const;

  /// @brief check and insert
  ///
  /// @param mState is the map to be filled
  /// @param surface is the surface to be checked for a Proxy
  void checkAndInsert(State& /*mState*/, const TrackingVolume& volume) const;

  /// @brief check and insert
  ///
  /// @param mState is the map to be filled
  /// @param surface is the surface to be checked for a Proxy
  void collectMaterialSurface(State& /*mState*/, const TrackingVolume& tVolume) const;

  /// Standard logger method
  const Logger& logger() const { return *m_logger; }

  /// The configuration object
  Config m_cfg;

  /// The straight line propagator
  StraightLinePropagator m_propagator;

  /// The logging instance
  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Acts

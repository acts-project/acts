// This file is part of the Acts project.
//
// Copyright (C) 2018-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/AccumulatedSurfaceMaterial.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Material/interface/IMaterialMapper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <map>
#include <memory>
#include <vector>

namespace Acts {

class ISurfaceMaterial;
class TrackingGeometry;
class TrackingVolume;

namespace Experimental {
class DetectorVolume;
class Detector;
}  // namespace Experimental

/// @brief selector for finding surface
struct MaterialSurface {
  bool operator()(const Surface& sf) const {
    return (sf.surfaceMaterial() != nullptr);
  }
};

/// @brief SurfaceMaterialMapper
///
/// This is the main feature tool to map material information
/// from a 3D geometry onto a surfaces based description.
///
/// It works for both, the `TrackingVolume` based detector description,
/// as well as the `DetectorVolume` one.
///
/// The process runs as such:
///
///  1) The geometry is parsed and for each Surface with
///     ProtoSurfaceMaterial a local store is initialized
///     the identification is done hereby through the
///     Surface::GeometryIdentifier
///
///  2) A Cache is generated that is used to keep the filling thread local,
///     the filling is protected with std::mutex
///
///  3) A number of N material tracks is read in, each track has :
///       origin, direction, material steps < position, step length, x0, l0, a,
///       z, rho >
///
///       For each track:
///          surfaces along the origin/direction path are collected
///          the closest material steps are assigned
///
///  4) Each 'hit' bin per event is counted and averaged at the end of the run
///
class SurfaceMaterialMapper final : public IMaterialMapper {
 public:
  using StraightLineTGPropagator = Propagator<StraightLineStepper, Navigator>;
  using StraightLineDetPropagator = 
    Propagator<StraightLineStepper, Experimental::DetectorNavigator>;

  /// @struct Config
  ///
  /// Nested Configuration struct for the material mapper
  struct Config {
    /// Correct for empty bins (recommended)
    bool emptyBinCorrection = true;
    /// Mapping output to debug stream
    bool mapperDebugOutput = false;
    /// Compute the variance of each material slab (only if using an input map)
    bool computeVariance = false;
    /// A potential veto
    MaterialInteractionVeto veto = NoVeto{};
  };

  /// Nested State struct which is used for the mapping prococess
  /// for caching the result
  ///
  /// In an appropriate way, the base class returned by the
  /// `createState(...)` interface method can be  transformed via
  /// `static_cast` tot he derived class needed by the mapper
  struct State final : public MaterialMappingState {
    /// @param [in] gctx The geometry context to use
    /// @param [in] mctx The magnetic field context to use
    State(const GeometryContext& gctx, const MagneticFieldContext& mctx)
        : geoContext(gctx), magFieldContext(mctx) {}

    /// The accumulated material per geometry ID
    std::map<GeometryIdentifier, AccumulatedSurfaceMaterial>
        accumulatedMaterial;

    /// The surface material of the input geometry
    std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>
        inputSurfaceMaterial;

    /// Reference to the geometry context for the mapping
    std::reference_wrapper<const GeometryContext> geoContext;

    /// Reference to the magnetic field context
    std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  };

  /// Delete the Default constructor
  SurfaceMaterialMapper() = delete;

  /// Constructor with config object
  ///
  /// @param cfg Configuration struct
  /// @param propagator The straight line propagator with the TrackingGeometry navigation
  /// @param slogger The logger
  SurfaceMaterialMapper(const Config& cfg,
                        StraightLineTGPropagator& propagator,
                        std::unique_ptr<const Logger> slogger =
                            getDefaultLogger("SurfaceMaterialMapper",
                                             Logging::INFO));

  /// Constructor with config object
  ///
  /// @param cfg Configuration struct
  /// @param propagator The straight line propagator with the Detector navigation
  /// @param slogger The logger
  SurfaceMaterialMapper(const Config& cfg,
                        StraightLineDetPropagator& propagator,
                        std::unique_ptr<const Logger> slogger =
                            getDefaultLogger("SurfaceMaterialMapper",
                                             Logging::INFO));

  /// @brief helper method that creates the cache for the mapping
  /// from a TrackingGeometry
  ///
  /// @param [in] gctx The geometry context to use
  /// @param [in] mctx The magnetic field context to use
  /// @param[in] tGeometry The geometry which should be mapped
  ///
  /// This method takes a TrackingGeometry,
  /// finds all surfaces with material proxis
  /// and returns you a Cache object to be used in the mapping process
  std::unique_ptr<MaterialMappingState> createState(
      const GeometryContext& gctx, const MagneticFieldContext& mctx,
      const TrackingGeometry& tGeometry) const final;

  /// @brief helper method that creates the cache for the mapping
  /// from a Detector object
  ///
  /// @param [in] gctx The geometry context to use
  /// @param [in] mctx The magnetic field context to use
  /// @param[in] detector The geometry which should be mapped
  ///
  /// This method takes a Detector object,
  /// finds all surfaces with material proxis
  /// and returns you a Cache object to be used in the mapping process
  std::unique_ptr<MaterialMappingState> createState(
      const GeometryContext& gctx, const MagneticFieldContext& mctx,
      const Experimental::Detector& detector) const final;

  /// Process/map a single track
  ///
  /// @param mState The current state map
  /// @param mTrack The material track to be mapped
  ///
  /// @note the RecordedMaterialSlab of the track are assumed
  /// to be ordered from the starting position along the starting direction
  ///
  /// The @return object is a pair of recorded material tracks, the mapped ones,
  /// and the unprocessed ones
  std::array<RecordedMaterialTrack, 2u> mapMaterialTrack(
      MaterialMappingState& mState,
      const RecordedMaterialTrack& mTrack) const final;

  /// @brief Method to finalize the maps
  ///
  /// It calls the final run averaging and then transforms
  /// the AccumulatedSurface material class to a surface material
  /// class type
  ///
  /// @param mState
  MaterialMappingResult finalizeMaps(MaterialMappingState& mState) const final;

 private:
  /// Loop through all the material interactions and add them to the
  /// associated surface
  ///
  /// @param mState The current state map
  /// @param mTrack The material track to be mapped
  ///
  std::array<RecordedMaterialTrack, 2u> mapInteraction(
      State& mState, const RecordedMaterialTrack& mTrack) const;

  /// Loop through all the material interactions and add them to the
  /// associated surface
  ///
  /// @note this option can be run to re-map material without the need of
  /// re-running the intersection, so it is primarily for optimisation in order
  /// to fine tune surface binning and description
  ///
  /// @param mState The current state map
  /// @param mTrack The material track to be mapped
  ///
  /// @note The material interactions are assumed to have an associated surface ID
  std::array<RecordedMaterialTrack, 2u> mapSurfaceInteraction(
      State& mState, const RecordedMaterialTrack& mTrack) const;

  /// @brief finds all surfaces with ProtoSurfaceMaterial of a tracking volume
  ///
  /// @param mState The state to be filled
  /// @param tVolume is current TrackingVolume
  void resolveMaterialSurfaces(State& mState,
                               const TrackingVolume& tVolume) const;

  /// @brief finds all surfaces with ProtoSurfaceMaterial of a detector volume
  ///
  /// @param mState The state to be filled
  /// @param dVolume is current DetectorVolume
  void resolveMaterialSurfaces(
      State& mState, const Experimental::DetectorVolume& dVolume) const;

  /// @brief Check and insert
  ///
  /// @param mState is the map to be filled
  /// @param surface is the surface to be checked for a Proxy
  void checkAndInsert(State& mState, const Surface& surface) const;

  /// Standard logger method
  const Logger& logger() const { return *m_logger; }

  /// The configuration object
  Config m_cfg;

  /// The straight line propagator
  std::shared_ptr<StraightLineTGPropagator> m_tgPropagator = nullptr;
  std::shared_ptr<StraightLineDetPropagator> m_detPropagator = nullptr;

  /// The logging instance
  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Acts

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/IAssignmentFinder.hpp"
#include "Acts/Material/ISurfaceMaterialAccumulator.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Material/MaterialInteractionAssignment.hpp"
#include "Acts/Material/TrackingGeometryMaterial.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <memory>
#include <span>
#include <utility>
#include <vector>

namespace Acts {

namespace detail {
class WorkerPool;
}  // namespace detail
/// Class that implements the material mapping procedure
/// @ingroup material_mapping
class MaterialMapper {
 public:
  /// @brief nested configuration struct
  struct Config {
    /// The assignment finder for material interaction assignments
    std::shared_ptr<const IAssignmentFinder> assignmentFinder = nullptr;
    /// The material accumulator for surfaces
    std::shared_ptr<const ISurfaceMaterialAccumulator>
        surfaceMaterialAccumulator = nullptr;
    /// Number of threads used by the batch @ref mapMaterial call, including
    /// the calling thread. The assignment finder must be safe to call
    /// concurrently when this is larger than one.
    std::size_t numThreads = 1;
  };

  /// @brief nested state struct
  ///
  /// It holds the states of the sub structs
  struct State {
    /// State of the surface material accumulator
    std::unique_ptr<ISurfaceMaterialAccumulator::State>
        surfaceMaterialAccumulatorState;
  };

  /// @brief nested options struct
  /// holds some options for the delegated calls
  struct Options {
    /// The assignment options (including vetos and re-assignments)
    MaterialInteractionAssignment::Options assignmentOptions;
  };

  /// @brief MaterialMapper constructor
  ///
  /// @param cfg the configuration struct
  /// @param mlogger the logger instance
  explicit MaterialMapper(
      const Config& cfg,
      std::unique_ptr<const Logger> mlogger =
          getDefaultLogger("BinnedSurfaceMaterialAccumulator", Logging::INFO));

  /// @brief Factory for creating the state
  /// @param gctx the geometry context
  /// @return Unique pointer to a new material mapping state object
  std::unique_ptr<State> createState(const GeometryContext& gctx) const;

  /// @brief Map the material interactions to the surfaces
  ///
  /// @param state the state object holding the sub states
  /// @param gctx the geometry context
  /// @param mctx the magnetic field context
  /// @param rmTrack the recorded material track
  /// @param options the call options (see above)
  ///
  /// @return the mapped and unmapped material tracks
  std::pair<RecordedMaterialTrack, RecordedMaterialTrack> mapMaterial(
      State& state, const GeometryContext& gctx,
      const MagneticFieldContext& mctx, const RecordedMaterialTrack& rmTrack,
      const Options& options = Options{}) const;

  /// @brief Map the material interactions of a batch of tracks
  ///
  /// The assignment of the tracks runs concurrently on the configured number
  /// of threads, the accumulation is delegated to the batch accumulation of
  /// the surface material accumulator. The result is identical to calling
  /// the single-track @ref mapMaterial for every track in order.
  ///
  /// @param state the state object holding the sub states
  /// @param gctx the geometry context
  /// @param mctx the magnetic field context
  /// @param rmTracks the recorded material tracks, which are not copied
  /// @param options the call options (see above)
  ///
  /// @return the mapped and unmapped material tracks, in input order
  std::vector<std::pair<RecordedMaterialTrack, RecordedMaterialTrack>>
  mapMaterial(State& state, const GeometryContext& gctx,
              const MagneticFieldContext& mctx,
              std::span<const RecordedMaterialTrack* const> rmTracks,
              const Options& options = Options{}) const;

  /// Finalize the maps
  /// @param state Material mapping state containing collected data
  /// @param gctx Geometry context for finalization
  /// @return Tracking geometry material map with finalized surface and volume materials
  TrackingGeometryMaterial finalizeMaps(const State& state,
                                        const GeometryContext& gctx) const;

 private:
  /// Access method to the logger
  const Logger& logger() const { return *m_logger; }

  /// The configuration
  Config m_cfg;

  /// The result of the assignment of a single track
  struct TrackAssignment {
    RecordedMaterialTrack mapped;
    RecordedMaterialTrack unmapped;
    std::vector<IAssignmentFinder::SurfaceAssignment> surfacesWithoutAssignment;
  };

  /// Assign the material of a single track, without accumulation
  TrackAssignment assignTrack(const GeometryContext& gctx,
                              const MagneticFieldContext& mctx,
                              const RecordedMaterialTrack& rmTrack,
                              const Options& options) const;

  /// The logger
  std::unique_ptr<const Logger> m_logger;

  /// The worker pool for batch mapping, nullptr for serial execution
  std::shared_ptr<detail::WorkerPool> m_pool;
};

}  // namespace Acts

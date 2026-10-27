// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Material/AccumulatedSurfaceMaterial.hpp"
#include "Acts/Material/ISurfaceMaterialAccumulator.hpp"
#include "Acts/Material/detail/MaterialSurfaceRegistry.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstdint>
#include <optional>
#include <unordered_map>

namespace Acts {

/// @brief The binned surface material accumulator
///
/// It consumes the assigned material interactions and then accumulates
/// the material on the surfaces in prepared binned containers for averaging
///
/// @ingroup material_mapping
class BinnedSurfaceMaterialAccumulator final
    : public ISurfaceMaterialAccumulator {
 public:
  /// @brief Nested config struct
  struct Config {
    /// Correct for empty bins (recommended)
    bool emptyBinCorrection = true;

    /// The surfaces to be used for the accumulation
    std::vector<const Surface*> materialSurfaces = {};
  };

  /// @brief Nested state struct
  struct State final : public ISurfaceMaterialAccumulator::State {
    /// The accumulated material per geometry ID
    std::map<GeometryIdentifier, AccumulatedSurfaceMaterial>
        accumulatedMaterial;

    /// Validated surface identities captured by createState.
    std::optional<detail::MaterialSurfaceRegistry> materialSurfaceRegistry;

    /// Cache for the batch accumulation: position of each geometry ID in
    /// @ref accumulatedMaterial, rebuilt when the number of entries changes
    std::unordered_map<GeometryIdentifier, std::uint32_t> slotIndex;
  };

  /// Constructor
  ///
  /// @param cfg the configuration struct
  /// @param mlogger the logger
  explicit BinnedSurfaceMaterialAccumulator(
      const Config& cfg,
      std::unique_ptr<const Logger> mlogger =
          getDefaultLogger("BinnedSurfaceMaterialAccumulator", Logging::INFO));

  /// Factory for creating the state
  /// @param gctx is the geometry context
  /// @return Unique pointer to newly created accumulator state
  std::unique_ptr<ISurfaceMaterialAccumulator::State> createState(
      const GeometryContext& gctx) const override;

  /// @brief Accumulate the material interaction on the surface
  ///
  /// Forwards to @ref accumulateTracks with a single track.
  /// @deprecated Use @ref accumulateTracks instead.
  ///
  /// @param state is the state of the accumulator
  /// @param gctx the geometry context
  /// @param interactions is the material interactions, with assigned surfaces
  /// @param surfacesWithoutAssignment are the surfaces without assignment
  [[deprecated("Use BinnedSurfaceMaterialAccumulator::accumulateTracks")]]
  void accumulate(ISurfaceMaterialAccumulator::State& state,
                  const GeometryContext& gctx,
                  const std::vector<MaterialInteraction>& interactions,
                  const std::vector<IAssignmentFinder::SurfaceAssignment>&
                      surfacesWithoutAssignment) const override;

  /// @brief Accumulate the material of a batch of tracks
  ///
  /// The accumulation is sharded by surface: every surface processes its
  /// share of the batch in track order, so the result does not depend on
  /// the worker pool, and the surfaces can be processed concurrently.
  ///
  /// @param state is the state of the accumulator
  /// @param gctx the geometry context
  /// @param tracks the per-track accumulation inputs, in mapping order
  /// @param pool optional worker pool, nullptr means serial execution
  void accumulateTracks(ISurfaceMaterialAccumulator::State& state,
                        const GeometryContext& gctx,
                        std::span<const TrackAccumulation> tracks,
                        detail::WorkerPool* pool) const override;

  /// Finalize the surface material maps
  ///
  /// @param state the state of the accumulator
  /// @param gctx is the geometry context
  ///
  /// @note this does the run average over the (binned) material
  /// @return Map of surface materials indexed by geometry identifiers
  std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>
  finalizeMaterial(ISurfaceMaterialAccumulator::State& state,
                   const GeometryContext& gctx) const override;

  /// Finalize assignments with the stable keys captured when creating state.
  /// @param state State created by this accumulator
  /// @param gctx Geometry context
  /// @return Finalized keyed and ID-only material maps
  TrackingGeometryMaterial finalizeMaps(
      ISurfaceMaterialAccumulator::State& state,
      const GeometryContext& gctx) const override;

 private:
  /// Access method to the logger
  const Logger& logger() const { return *m_logger; }

  /// The configuration
  Config m_cfg;

  /// The logger
  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Acts

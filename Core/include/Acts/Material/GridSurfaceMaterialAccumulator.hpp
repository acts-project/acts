// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Material/AccumulatedMaterialSlab.hpp"
#include "Acts/Material/interface/ISurfaceMaterialAccumulator.hpp"
#include "Acts/Utilities/IMultiAxis.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <map>
#include <memory>
#include <vector>

namespace Acts {

/// @brief The grid surface material accumulator
///
/// It consumes the assigned material interactions and then accumulates
/// the material on the surfaces in a 2D grid (see @c GridSurfaceMaterial),
/// exactly as @c BinnedSurfaceMaterialAccumulator does for @c BinUtility -
/// this is its modern counterpart, producing @c GridSurfaceMaterial instead
/// of @c BinnedSurfaceMaterial.
///
/// Surfaces are expected to carry either a @c ProtoGridSurfaceMaterial
/// (whose possibly deferred axes are resolved against the surface, exactly
/// as @c BinnedSurfaceMaterialAccumulator resolves one against a
/// @c BinUtility) or an already concrete @c GridSurfaceMaterial to be
/// re-accumulated; surfaces carrying neither fall back to homogeneous
/// accumulation, producing a @c HomogeneousSurfaceMaterial. The legacy
/// @c BinUtility-based @c ProtoSurfaceMaterial and @c BinnedSurfaceMaterial
/// are not consumed here - that remains
/// @c BinnedSurfaceMaterialAccumulator's job.
class GridSurfaceMaterialAccumulator final
    : public ISurfaceMaterialAccumulator {
 public:
  /// The storage backend the finalized @c GridSurfaceMaterial is built with.
  ///
  /// @c GloballyIndexed storage requires knowing which bins, across which
  /// surfaces, should share a material vector entry - that is a
  /// postprocessing step run over the finalized maps, not something this
  /// accumulator can decide on its own, so it is not offered here.
  enum class StorageKind { Direct, Indexed };

  /// @brief Nested config struct
  struct Config {
    /// Correct for empty bins (recommended)
    bool emptyBinCorrection = true;

    /// The surfaces to be used for the accumulation
    std::vector<const Surface*> materialSurfaces = {};

    /// The storage backend the finalized GridSurfaceMaterial is built with
    StorageKind storageKind = StorageKind::Indexed;
  };

  /// @brief Nested state struct
  struct State final : public ISurfaceMaterialAccumulator::State {
    /// Per-surface accumulation on a resolved 2D grid
    struct GridAccumulation {
      /// The multi-axis doing local position -> index, resolved once at
      /// createState() time
      std::unique_ptr<IMultiAxis2D> multiAxis;
      /// One accumulated slab per bin, including under-/overflow, flattened
      /// in multiAxis's global bin order
      std::vector<AccumulatedMaterialSlab> accumulatedMaterial;
    };

    /// The accumulated material per geometry ID, for surfaces with a
    /// resolved 2D grid
    std::map<GeometryIdentifier, GridAccumulation> gridMaterial;

    /// The accumulated material per geometry ID, for surfaces without a
    /// grid (homogeneous material)
    std::map<GeometryIdentifier, AccumulatedMaterialSlab> homogeneousMaterial;
  };

  /// Constructor
  ///
  /// @param cfg the configuration struct
  /// @param mlogger the logger
  explicit GridSurfaceMaterialAccumulator(
      const Config& cfg,
      std::unique_ptr<const Logger> mlogger =
          getDefaultLogger("GridSurfaceMaterialAccumulator", Logging::INFO));

  /// Factory for creating the state
  /// @param gctx is the geometry context
  /// @return Unique pointer to newly created accumulator state
  /// @throws std::invalid_argument if a configured surface has no material
  ///         assigned, or carries a legacy BinUtility-based
  ///         ProtoSurfaceMaterial/BinnedSurfaceMaterial
  std::unique_ptr<ISurfaceMaterialAccumulator::State> createState(
      const GeometryContext& gctx) const override;

  /// @brief Accumulate the material interaction on the surface
  ///
  /// @param state is the state of the accumulator
  /// @param gctx is the geometry context
  /// @param interactions is the material interactions, with assigned surfaces
  /// @param surfacesWithoutAssignment are the surfaces without assignment
  ///
  /// @note this the track average over the binned material
  void accumulate(ISurfaceMaterialAccumulator::State& state,
                  const GeometryContext& gctx,
                  const std::vector<MaterialInteraction>& interactions,
                  const std::vector<IAssignmentFinder::SurfaceAssignment>&
                      surfacesWithoutAssignment) const override;

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

 private:
  /// Resolve the local (loc0, loc1) position of a global position on a
  /// surface, in the surface's canonical local frame (i.e. the same frame
  /// resolveMultiAxis resolves a ProtoGridSurfaceMaterial's axes against)
  ///
  /// @param gctx the geometry context
  /// @param surface the surface
  /// @param position the global position
  /// @param direction the global direction
  /// @return the local 2D position
  /// @throws std::invalid_argument if the position cannot be resolved onto
  ///         the surface
  static Vector2 resolveLocalPosition(const GeometryContext& gctx,
                                      const Surface& surface,
                                      const Vector3& position,
                                      const Vector3& direction);

  /// Access method to the logger
  const Logger& logger() const { return *m_logger; }

  /// The configuration
  Config m_cfg;

  /// The logger
  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Acts

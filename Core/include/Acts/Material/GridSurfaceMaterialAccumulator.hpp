// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/AccumulatedMaterialSlab.hpp"
#include "Acts/Material/detail/MaterialSurfaceRegistry.hpp"
#include "Acts/Material/interface/ISurfaceMaterialAccumulator.hpp"
#include "Acts/Utilities/IMultiAxis.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <map>
#include <memory>
#include <optional>
#include <vector>

namespace Acts {

/// Accumulate assigned interactions into two-dimensional grid surface material.
///
/// Prototype axes are resolved against the surface. Existing grid axes are
/// reused for remapping; homogeneous material uses two one-bin axes. Global
/// intersections are converted to surface-local coordinates using the geometry
/// context. Each touched bin is averaged once per track, then over all tracks.
/// Finalized grids use direct storage, including underflow and overflow bins,
/// and the same split factor (zero) as BinnedSurfaceMaterialAccumulator.
///
/// @ingroup material_mapping
class GridSurfaceMaterialAccumulator final
    : public ISurfaceMaterialAccumulator {
 public:
  /// Accumulator configuration.
  struct Config {
    /// Include crossings without material in the per-bin track average.
    bool emptyBinCorrection = true;
    /// Surfaces carrying prototype, grid or homogeneous material.
    std::vector<const Surface*> materialSurfaces;
  };

  /// Per-run accumulation state.
  struct State final : public ISurfaceMaterialAccumulator::State {
    /// Accumulation data for one surface.
    struct AccumulatedGrid {
      /// Surface used to convert global intersections to local coordinates.
      const Surface* surface = nullptr;
      /// Resolved axes in material storage order.
      std::unique_ptr<IMultiAxis2D> multiAxis;
      /// Per-bin accumulation data, including underflow and overflow bins.
      std::vector<AccumulatedMaterialSlab> material;
      /// Whether the material axes are reversed relative to the surface axes.
      bool swapAxes = false;
    };

    /// Accumulated grids indexed by geometry identifier.
    std::map<GeometryIdentifier, AccumulatedGrid> accumulatedMaterial;
    /// Validated identities captured when creating this state.
    std::optional<detail::MaterialSurfaceRegistry> materialSurfaceRegistry;
  };

  /// Construct an accumulator.
  /// @param cfg Accumulator configuration
  /// @param mlogger Logger
  explicit GridSurfaceMaterialAccumulator(
      const Config& cfg,
      std::unique_ptr<const Logger> mlogger =
          getDefaultLogger("GridSurfaceMaterialAccumulator", Logging::INFO));

  /// @copydoc ISurfaceMaterialAccumulator::createState
  std::unique_ptr<ISurfaceMaterialAccumulator::State> createState(
      const GeometryContext& gctx) const override;

  /// @copydoc ISurfaceMaterialAccumulator::accumulate
  void accumulate(ISurfaceMaterialAccumulator::State& state,
                  const GeometryContext& gctx,
                  const std::vector<MaterialInteraction>& interactions,
                  const std::vector<IAssignmentFinder::SurfaceAssignment>&
                      surfacesWithoutAssignment) const override;

  /// @copydoc ISurfaceMaterialAccumulator::finalizeMaterial
  std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>
  finalizeMaterial(ISurfaceMaterialAccumulator::State& state,
                   const GeometryContext& gctx) const override;

  /// @copydoc ISurfaceMaterialAccumulator::finalizeMaps
  TrackingGeometryMaterial finalizeMaps(
      ISurfaceMaterialAccumulator::State& state,
      const GeometryContext& gctx) const override;

 private:
  /// Access the logger.
  const Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Acts

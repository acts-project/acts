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
#include "Acts/Material/interface/ISurfaceMaterialAccumulator.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

/// @brief The binned surface material accumulator
///
/// It consumes the assigned material interactions and then accumulates
/// the material on the surfaces in prepared binned containers for averaging

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
  /// Access method to the logger
  const Logger& logger() const { return *m_logger; }

  /// The configuration
  Config m_cfg;

  /// The logger
  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Acts

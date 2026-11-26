// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Vertexing/IVertexFinder.hpp"
#include "Acts/Vertexing/IVertexFitter.hpp"

#include <memory>

namespace Acts {

/// Minimal pairwise vertex finder for displaced two-track candidates.
/// Fits each accepted pair independently, without overlap removal, mass
/// tagging, primary-vertex selection or a displacement cut. Complexity is
/// quadratic in the number of input tracks. Intended as a seed-and-fit building
/// block. Seeding approximates the field by its axial component at the pair
/// midpoint; non-axial fields are rejected. Fitting uses the supplied fitter's
/// propagation.
class TwoTrackVertexFinder final : public IVertexFinder {
 public:
  /// Finder configuration.
  struct Config {
    /// Single-vertex fitter, required.
    std::shared_ptr<const IVertexFitter> vertexFitter;
    /// Field provider for the reusable query, required.
    std::shared_ptr<const MagneticFieldProvider> field;
    /// Extract parameters from input tracks, required.
    InputTrack::Extractor extractParameters;
    /// Reject seeds with larger trajectory-point separation.
    double maxSeedSeparation = 1. * UnitConstants::mm;
    /// Reject fits with larger total chi-squared.
    double maxChi2 = 25.;
    /// Restrict pairs to opposite charges.
    bool oppositeChargeOnly = true;
  };

  /// @param cfg Configuration, whose tools must support the same track type
  /// @throws std::invalid_argument for missing tools or invalid cuts
  explicit TwoTrackVertexFinder(Config cfg);

  /// Fit independent two-track candidates.
  /// @param tracks Input tracks, which must outlive the returned vertices
  /// @param options Contexts; useConstraintInFit must be false
  /// @param state State from makeState for options.magFieldContext
  /// @return Accepted vertices or a field, input or fitting error
  /// Seed degeneracies reject the pair; field and fitter errors are propagated.
  /// Fit constraints are rejected: this finder performs unconstrained fits.
  Result<std::vector<Vertex>> find(const std::vector<InputTrack>& tracks,
                                   const VertexingOptions& options,
                                   IVertexFinder::State& state) const override;

  /// Create a reusable query and track-removal state.
  /// @param mctx Borrowed field context, which must outlive the state
  /// @return State for this finder
  IVertexFinder::State makeState(
      const MagneticFieldContext& mctx) const override;

  /// @param state State from makeState
  /// @param tracks Tracks to exclude
  /// Excludes these tracks on subsequent calls using this state.
  void setTracksToRemove(IVertexFinder::State& state,
                         const std::vector<InputTrack>& tracks) const override;

 private:
  Config m_cfg;
};
}  // namespace Acts

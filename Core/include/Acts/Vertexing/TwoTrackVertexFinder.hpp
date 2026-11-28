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
#include "Acts/Vertexing/TrackAtVertex.hpp"

#include <memory>

namespace Acts {

class IVertexFitter;

/// @ingroup vertexing
/// Minimal pairwise vertex finder for displaced two-track candidates.
/// Fits each accepted pair independently through IVertexFitter, with an
/// optional per-call position constraint. This is a reconstruction building
/// block, not a production secondary-vertex selection.
///
/// @par Usage
/// @code{.cpp}
/// Acts::TwoTrackVertexFinder::Config cfg;
/// cfg.vertexFitter = fitter;  // shared_ptr<const IVertexFitter>
/// cfg.field = fieldProvider;  // shared_ptr<const MagneticFieldProvider>
/// cfg.extractParameters.connect<&Acts::InputTrack::extractParameters>();
/// Acts::TwoTrackVertexFinder finder(cfg);
/// auto state = finder.makeState(mctx);
/// Acts::VertexingOptions options(gctx, mctx);  // Unconstrained fit.
/// auto result = finder.find(inputTracks, options, state);
/// @endcode
///
/// @par Seeding and selection
/// The internal estimator uses transverse circle intersections and longitudinal
/// track separation to choose a seed. For nonintersecting circles it uses their
/// closest transverse points; in zero field it uses closest points on straight
/// lines. Its point separation is not an exact three-dimensional helix DCA.
/// Only the nearest half-turn of each trajectory is considered. Parallel lines
/// and concentric circles without a unique seed reject the pair.
///
/// The field is sampled at the midpoint of the input track positions and
/// approximated as uniform and axial for seeding. Non-axial fields are
/// rejected; the fitter remains responsible for propagation. Its propagator and
/// the finder's field provider must describe the same field.
///
/// The default selection requires opposite charges, seed separation at most
/// 1 mm, and total fit chi-squared at most 25. These are prototype cuts, not
/// calibrated physics working points. There is no displacement cut, mass
/// tagging, overlap resolution or vertex merging. Output vertices may share
/// tracks. An explicitly enabled VertexingOptions constraint is passed to every
/// pair's fit without changing its geometric seed or seed-separation
/// preselection. The total fit chi-squared cut then includes the constraint
/// contribution. If disabled, the constraint is ignored, including its position
/// and covariance. A constraint requires a finite position and a finite,
/// symmetric, positive definite full 4D covariance; a spatial-only covariance
/// with zero time variance is rejected. It is a Gaussian position/time prior,
/// not a pointing constraint or a reference PV for selection. A prior
/// reconstructed from the input tracks is not an independent measurement; this
/// API does not model such correlations. Runtime scales quadratically with the
/// input multiplicity.
///
/// @par State and failures
/// Keep the field context alive and unchanged while using the state, and use
/// the same context in VertexingOptions. Construct a fresh state for each
/// event. setTracksToRemove() excludes tracks on subsequent calls with that
/// state. Inputs are deduplicated by track identity. Seed degeneracies reject
/// the pair; other seed errors, field errors and fitter errors fail the call.
///
/// @par Synthetic validation
/// The TwoTrackVertexFinder unit-test executable generates correlated pion
/// pairs from a 497.611 MeV parent, boosted to momenta of 1, 2 and 5 GeV, at
/// transverse radii of 1, 10, 50 and 100 mm. The parent moves transversely;
/// vertex z is uniform within +/-20 mm, and the rest-frame angular sample has
/// |cos(theta)| < 0.8. There are 20 pairs per radius/momentum bin.
///
/// An Eigen stepper propagates the daughters through a uniform axial 2 T field
/// to a common origin perigee. Unsmeared seeds must recover the known
/// production point. Independent Gaussian smearing then assumes 20 micrometers
/// in d0, 50 micrometers in z0, 0.2 mrad in each angle, 0.1 percent in q/p and
/// 0.1 ns in time. Billoir receives these diagonal covariances, but neither the
/// true vertex nor a mass constraint.
///
/// Run the scan and print its per-bin acceptance and aggregate residual, pull
/// and mass diagnostics with:
/// @code{.sh}
/// cmake --build build --target ActsUnitTestTwoTrackVertexFinder
/// build/bin/ActsUnitTestTwoTrackVertexFinder --log_level=message
/// @endcode
/// The fixed-sample regression checks require at least 16/20 accepted pairs in
/// each bin, spatial residual RMS below 2 mm, mean mass within 3 MeV of the
/// injected mass, and spatial pull RMS between 0.5 and 1.5. Other tests cover
/// analytic geometry, interface dispatch, error forwarding and track exclusion.
///
/// This track-level sample has no parent transport, physical decay-radius
/// distribution, detector hits, material, pattern recognition or realistic
/// combinatorial background. Acceptance measures seed-and-fit behavior for
/// supplied tracks, not detector efficiency or a calibration of uncertainties.
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
  /// @param options Contexts and optional per-call position constraint
  /// @param state State from makeState for options.magFieldContext
  /// @return Accepted vertices or a field, input or fitting error
  /// Seed degeneracies reject the pair; field and fitter errors are propagated.
  /// Malformed enabled constraints are rejected before fitting any pairs.
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

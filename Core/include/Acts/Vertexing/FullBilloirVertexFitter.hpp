// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/IVertexFitter.hpp"
#include "Acts/Vertexing/TrackLinearizer.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"

#include <span>

namespace Acts {

/// @class FullBilloirVertexFitter
///
/// @brief Vertex fitter class implementing the Billoir vertex fitter
///
/// This class implements the Billoir vertex fitter from Ref. (1). It is also
/// useful to have a look at Ref. (2). The cross-covariance matrices are derived
/// in Ref. (3). Note that the Billoir vertex fitter outputs one 4D vertex
/// position and nTrack momenta at this very point.
///
/// Ref. (1):
/// Fast vertex fitting with a local parametrization of tracks.
/// Author(s) Billoir, P ; Qian, S
/// In: Nucl. Instrum. Methods Phys. Res., A 311 (1992) 139-150
/// DOI 10.1016/0168-9002(92)90859-3
///
/// Ref. (2):
/// Pattern Recognition, Tracking and Vertex Reconstruction in Particle
/// Detectors.
/// Author(s) Fruehwirth, R ; Strandli, A
///
/// Ref. (3):
/// ACTS White Paper: Cross-Covariance Matrices in the Billoir Vertex Fit
/// https://acts.readthedocs.io/en/latest/white_papers/billoir-covariances.html
/// Author(s) Russo, F
class FullBilloirVertexFitter final : public IVertexFitter {
 public:
  /// Configuration options for the Billoir vertex fitter.
  struct Config {
    /// Maximum number of iterations in fitter
    int maxIterations = 5;

    /// Function to extract parameters from InputTrack
    InputTrack::Extractor extractParameters;

    /// Track linearizer
    TrackLinearizer trackLinearizer;

    /// Magnetic field provider, used only to create the fitter cache in
    /// @c makeCache. This has to be the same field that the track linearizer
    /// uses.
    ///
    /// Optional: it is required only to drive this fitter through the
    /// @c IVertexFitter interface. Callers using @c fit directly supply their
    /// own field cache and can leave this unset.
    std::shared_ptr<const MagneticFieldProvider> bField;
  };

  /// @brief The fitter-private cache
  ///
  /// The Billoir fit keeps no state across vertices; the only scratch data is
  /// the magnetic field cache used when linearizing tracks.
  struct Cache {
    /// Constructor for the Billoir fitter cache
    /// @param field Magnetic field provider for track extrapolation
    /// @param magContext Magnetic field context for field evaluations
    Cache(const MagneticFieldProvider& field,
          const Acts::MagneticFieldContext& magContext);

    /// Magnetic field cache for field evaluations during fitting
    MagneticFieldProvider::Cache fieldCache;
  };

  /// @brief Constructor for user-defined InputTrack type
  ///
  /// @param cfg Configuration object
  /// @param logger Logging instance
  explicit FullBilloirVertexFitter(
      const Config& cfg,
      std::unique_ptr<const Logger> logger =
          getDefaultLogger("FullBilloirVertexFitter", Logging::INFO));

  /// @brief Fit method, fitting vertex for provided tracks with constraint
  ///
  /// @param paramVector Vector of track objects to fit vertex to
  /// @param vertexingOptions Vertexing options
  /// @param fieldCache The magnetic field cache
  ///
  /// @return Fitted vertex
  Result<Vertex> fit(std::span<const InputTrack> paramVector,
                     const VertexingOptions& vertexingOptions,
                     MagneticFieldProvider::Cache& fieldCache) const;

  /// @copydoc IVertexFitter::makeCache
  ///
  /// @note Requires @c Config::bField to be set. It is not needed when calling
  /// @c fit directly with a caller-supplied field cache, which is why it is
  /// checked here rather than in the constructor.
  IVertexFitter::Cache makeCache(
      const MagneticFieldContext& mctx) const override;

  /// @copydoc IVertexFitter::fit
  ///
  /// Each vertex in @p problem is fitted independently: unlike the adaptive
  /// multi-vertex fitter, tracks are not competed between vertices.
  Result<void> fit(VertexFitProblem& problem,
                   const VertexingOptions& vertexingOptions,
                   IVertexFitter::Cache& cache) const override;

  /// @copydoc IVertexFitter::addVertices
  ///
  /// Since vertices are fitted independently, only @p newVertices are fitted;
  /// vertices already in @p problem are left untouched.
  Result<void> addVertices(VertexFitProblem& problem,
                           std::span<Vertex* const> newVertices,
                           const VertexingOptions& vertexingOptions,
                           IVertexFitter::Cache& cache) const override;

  /// @copydoc IVertexFitter::fitSingle
  Result<Vertex> fitSingle(std::span<const InputTrack> trackVector,
                           const VertexingOptions& vertexingOptions,
                           IVertexFitter::Cache& cache) const override;

 private:
  /// @brief Fits a single vertex of the problem in place
  ///
  /// @param problem The multi-vertex fit problem
  /// @param vtx The vertex to fit
  /// @param vertexingOptions Vertexing options
  /// @param cache Fitter cache
  Result<void> fitVertex(VertexFitProblem& problem, Vertex* vtx,
                         const VertexingOptions& vertexingOptions,
                         Cache& cache) const;

  /// Configuration object
  Config m_cfg;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to logging instance
  const Logger& logger() const;
};

}  // namespace Acts

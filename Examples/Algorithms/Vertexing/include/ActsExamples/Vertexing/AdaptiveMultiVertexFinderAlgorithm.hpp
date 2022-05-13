// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"

#include <string>

namespace ActsExamples {

class AdaptiveMultiVertexFinderAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// Input track parameters collection
    std::string inputTrackParameters;
    /// Output proto vertex collection
    std::string outputProtoVertices;
    /// Output vertex collection
    std::string outputVertices = "vertices";
    /// Output reconstruction time in ms
    std::string outputTime = "time";
    /// The magnetic field
    std::shared_ptr<Acts::MagneticFieldProvider> bField;

    //beam spot constraint
    bool amvf_useBeamSpotConstraint = true;

    // Max z interval used for adding tracks to fit:
    float amvf_tracksMaxZinterval = 3. * Acts::UnitConstants::mm;

    // Maximum allowed significance of track position to vertex seed
    // to consider track as compatible track for vertex fit
    float amvf_tracksMaxSignificance = 5.;

    // Max chi2 value for which tracks are considered compatible with
    // the fitted vertex. These tracks are removed from the seedTracks
    // after the fit has been performed.
    float amvf_maxVertexChi2 = 18.42;

    // Maximum significance on the distance between two vertices
    // to allow merging of two vertices.
    float amvf_maxMergeVertexSignificance = 3.;

    // Minimum weight a track has to have to be considered a compatible
    // track with a vertex candidate.
    // Note: This value has to be the same as the one in the AMVFitter.
    float amvf_minWeight = 0.0001;

    // Maximal number of iterations in the finding procedure
    int amvf_maxIterations = 100;

    // Maximum vertex contamination value
    float amvf_maximumVertexContamination = 0.5;

    // Use the full available vertex covariance information after
    // seeding for the IP estimation. In original implementation
    // this is not (!) done, however, this is probably not correct.
    // So definitely consider setting this to true.
    bool amvf_useVertexCovForIPEstimation = false;

  };

  AdaptiveMultiVertexFinderAlgorithm(const Config& config,
                                     Acts::Logging::Level level);

  /// Find vertices using the adapative multi vertex finder algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples

// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/TelescopeSeedingAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Seeding/EstimateTrackParamsFromSeed.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <cmath>
#include <cstddef>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <system_error>
#include <utility>
#include <vector>

ActsExamples::TelescopeSeedingAlgorithm::TelescopeSeedingAlgorithm(
    ActsExamples::TelescopeSeedingAlgorithm::Config cfg,
    Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("TelescopeSeedingAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurements input collection");
  }
  if (m_cfg.outputTrackParameters.empty()) {
    throw std::invalid_argument("Missing track parameters output collection");
  }
  if (!m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }

  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);

  // Set up the track parameters covariance (the same for all tracks)
  for (std::size_t i = Acts::eBoundLoc0; i < Acts::eBoundSize; ++i) {
    m_covariance(i, i) = m_cfg.initialVarInflation[i] * m_cfg.initialSigmas[i] *
                         m_cfg.initialSigmas[i];
  }
}

ActsExamples::ProcessCode ActsExamples::TelescopeSeedingAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  auto const& measurements = m_inputMeasurements(ctx);
  ACTS_VERBOSE("Read " << measurements.size() << " measurements");

  TrackParametersContainer trackParameters;

  IndexSourceLink::SurfaceAccessor surfaceAccessor{*m_cfg.trackingGeometry};

  // Loop over all measurements to estimate track parameters
  for (std::size_t imeas = 0; imeas < measurements.size(); ++imeas) {
    const auto& measurement = measurements[imeas];

    std::visit(
        [&](const auto& m) {
          const Acts::Surface* surface = surfaceAccessor(m.sourceLink());

          if (surface == nullptr) {
            ACTS_WARNING(
                "Surface from source link is not found in the tracking "
                "geometry");
            return true;
          }

          // The measurement on the first layer
          if (surface->geometryId().layer() != m_cfg.selectedLayer) {
            return true;
          }

          // The parameters and covariance for the measured position
          const auto& measPar = m.parameters();
          const auto& measCov = m.covariance();
          static constexpr auto measDim = m.size();

          Acts::BoundVector params = Acts::BoundVector::Zero();
          Acts::BoundSquareMatrix cov = m_covariance;
          params[Acts::eBoundTheta] = m_cfg.theta;
          params[Acts::eBoundPhi] = m_cfg.phi;
          params[Acts::eBoundQOverP] = m_cfg.qOp;

          if (m.contains(Acts::eBoundLoc0)) {
            params[Acts::eBoundLoc0] = measPar[Acts::eBoundLoc0];
            cov(Acts::eBoundLoc0, Acts::eBoundLoc0) =
                measCov(Acts::eBoundLoc0, Acts::eBoundLoc0);
          }
          if (m.contains(Acts::eBoundLoc1)) {
            params[Acts::eBoundLoc1] = measPar[Acts::eBoundLoc1];
            cov(Acts::eBoundLoc1, Acts::eBoundLoc1) =
                measCov(Acts::eBoundLoc1, Acts::eBoundLoc1);
          }
          if (m.contains(Acts::eBoundTime)) {
            params[Acts::eBoundTime] = measPar[measDim - 1];
            cov(Acts::eBoundTime, Acts::eBoundTime) =
                measCov(measDim - 1, measDim - 1);
          } else {
            params[Acts::eBoundTime] = m_cfg.time;
          }

          trackParameters.emplace_back(surface->getSharedPtr(), params, cov,
                                       m_cfg.particleHypothesis);
          return true; 
	},
        measurement);
  }

  ACTS_VERBOSE("Estimated " << trackParameters.size() << " track parameters");

  m_outputTrackParameters(ctx, std::move(trackParameters));

  return ProcessCode::SUCCESS;
}

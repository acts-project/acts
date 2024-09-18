// This file is part of the Acts project.
//
// Copyright (C) 2021-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
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

namespace ActsExamples {

namespace {

Acts::BoundSquareMatrix makeInitialCovariance(
    const TrackParamsEstimationAlgorithm::Config& config,
    const Acts::BoundVector& params, const SimSpacePoint& sp) {
  Acts::BoundSquareMatrix result = Acts::BoundSquareMatrix::Zero();

  for (std::size_t i = Acts::eBoundLoc0; i < Acts::eBoundSize; ++i) {
    double sigma = config.initialSigmas[i];

    // Add momentum dependent uncertainties
    sigma +=
        config.initialSimgaQoverPCoefficients[i] * params[Acts::eBoundQOverP];

    double var = sigma * sigma;

    // Inflate the time uncertainty if no time measurement is available
    if (i == Acts::eBoundTime && !sp.t().has_value()) {
      var *= config.noTimeVarInflation;
    }

    // Inflate the initial covariance
    var *= config.initialVarInflation[i];

    result(i, i) = var;
  }

  return result;
}

}  // namespace

ActsExamples::TrackParamsEstimationAlgorithm::TrackParamsEstimationAlgorithm(
    ActsExamples::TrackParamsEstimationAlgorithm::Config cfg,
    Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("TrackParamsEstimationAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputSeeds.empty()) {
    throw std::invalid_argument("Missing seeds input collection");
  }
  if (m_cfg.outputTrackParameters.empty()) {
    throw std::invalid_argument("Missing track parameters output collection");
  }
  if (!m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }
  if (!m_cfg.magneticField) {
    throw std::invalid_argument("Missing magnetic field");
  }

  m_inputSeeds.initialize(m_cfg.inputSeeds);
  m_inputTracks.maybeInitialize(m_cfg.inputProtoTracks);

  m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);
  m_outputSeeds.maybeInitialize(m_cfg.outputSeeds);
  m_outputTracks.maybeInitialize(m_cfg.outputProtoTracks);
}

ActsExamples::ProcessCode ActsExamples::TrackParamsEstimationAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  auto const& seeds = m_inputSeeds(ctx);

  //==================================================================
  // MODIFIED BY ME
  //      // rotation of SP should be done here
  //      const auto bottomSP = seed.sp().at(0);
  //      const auto middleSP = seed.sp().at(1);
  //      const auto topSP = seed.sp().at(2);
  //
  //      std::cout << "before rotation" << std::endl;
  //      std::cout << "SPinseed_bottomsp= " << bottomSP->x() << " " <<
  //      bottomSP->y() << " " << bottomSP->z() << " " << std::endl; std::cout
  //      << "SPinseed_middlesp= " << middleSP->x() << " " << middleSP->y() <<
  //      " " << middleSP->z() << " " << std::endl; std::cout <<
  //      "SPinseed_topsp= " << topSP->x() << " " << topSP->y() << " " <<
  //      topSP->z() << " " << std::endl;
  //
  //      double_t ytmpb = bottomSP->y();
  //      double_t ytmpm = middleSP->y();
  //      double_t ytmpt = topSP->y();
  //      bottomSP->y() == bottomSP->z();
  //      bottomSP->z() == -ytmpb;
  //      middleSP->y() == middleSP->z();
  //      middleSP->z() == -ytmpm;
  //      topSP->y() == topSP->z();
  //      topSP->z() == -ytmpt;
  //
  //      std::cout << "after rotation" << std::endl;
  //      std::cout << "bottomsp= " << bottomSP->x() << " " << bottomSP->y()
  //      << " " << bottomSP->z() << std::endl; std::cout << "middlesp= " <<
  //      middleSP->x() << " " << middleSP->y() << " " << middleSP->z() <<
  //      std::endl; std::cout << "topsp= " << topSP->x() << " " << topSP->y()
  //      << " " << topSP->z() <<  std::endl;
  //===================================================================

  ACTS_VERBOSE("Read " << seeds.size() << " seeds");

  TrackParametersContainer trackParameters;
  trackParameters.reserve(seeds.size());

  SimSeedContainer outputSeeds;
  if (m_outputSeeds.isInitialized()) {
    outputSeeds.reserve(seeds.size());
  }

  const ProtoTrackContainer* inputTracks = nullptr;
  ProtoTrackContainer outputTracks;
  if (m_inputTracks.isInitialized() && m_outputTracks.isInitialized()) {
    const auto& inputTracksRef = m_inputTracks(ctx);
    if (seeds.size() != inputTracksRef.size()) {
      ACTS_FATAL("Inconsistent number of seeds and prototracks");
      return ProcessCode::ABORT;
    }
    inputTracks = &inputTracksRef;
    outputTracks.reserve(seeds.size());
  }

  auto bCache = m_cfg.magneticField->makeCache(ctx.magFieldContext);

  IndexSourceLink::SurfaceAccessor surfaceAccessor{*m_cfg.trackingGeometry};

  // Loop over all found seeds to estimate track parameters
  for (std::size_t iseed = 0; iseed < seeds.size(); ++iseed) {
    const auto& seed = seeds[iseed];
    // Get the bottom space point and its reference surface
    const auto bottomSP = seed.sp().front();
    const auto middleSP = seed.sp()[1];
    if (bottomSP->sourceLinks().empty()) {
      ACTS_WARNING("Missing source link in the space point");
      continue;
    }
    const auto& sourceLink = bottomSP->sourceLinks()[0];
    const Acts::Surface* surface = surfaceAccessor(sourceLink);

    if (surface == nullptr) {
      ACTS_WARNING(
          "Surface from source link is not found in the tracking geometry");
      continue;
    }

    if (m_cfg.verbose)
      std::cout << "TrackParamsEstimationAlgorithm_bottomSP= " << bottomSP->x()
                << " " << bottomSP->y() << " " << bottomSP->z() << std::endl;


    // Get the magnetic field at the bottom space point
    auto fieldRes = m_cfg.magneticField->getField(
        {middleSP->x(), middleSP->z(), -middleSP->y()}, bCache);
    if (!fieldRes.ok()) {
      ACTS_ERROR("Field lookup error: " << fieldRes.error());
      return ProcessCode::ABORT;
    }
    Acts::Vector3 field = *fieldRes;

    if (!m_cfg.truthSeeding) {
      field.y() = -field.z();
      field.z() = 0.0;
    }

    // Estimate the track parameters from seed
    auto optParams = Acts::estimateTrackParamsFromSeed(
        ctx.geoContext, seed.sp().begin(), seed.sp().end(), *surface, field,
         m_cfg.bFieldMin, logger(), m_cfg.verbose);

    // ===================== ALTERNATIVE WAY OF EXTRACTING TRACK PARAMETERS
    // (Noemi)
    // - SP SHOULD BE ROTATED (where???)
    // - B FIELD SHOULD NOT BE ROTATED
    // - ALTERNATIVE estimateTrackParameter IS USED
    //
    //    std::cout << "TrackParamsEstimationAlgorithm_bottomSP= " <<
    //    bottomSP->x() << " " << bottomSP->y() << " " << bottomSP->z() <<
    //    std::endl; std::cout << "fieldB= " << field.x() << " " << field.y() <<
    //    " " << field.z() << std::endl; auto optParams =
    //    Acts::estimateTrackParamsFromSeed( seed.sp().begin(),
    //    seed.sp().end(),logger());
    //

    // END OF CHANGES
    //===============================================================================

    if (!optParams.has_value()) {
      ACTS_WARNING("Estimation of track parameters for seed " << iseed
                                                              << " failed.");
      continue;
    }

    const auto& params = optParams.value();

    Acts::BoundSquareMatrix cov =
        makeInitialCovariance(m_cfg, params, *bottomSP);

    trackParameters.emplace_back(surface->getSharedPtr(), params, cov,
                                 m_cfg.particleHypothesis);
    if (m_outputSeeds.isInitialized()) {
      outputSeeds.push_back(seed);
    }
    if (m_outputTracks.isInitialized() && inputTracks != nullptr) {
      outputTracks.push_back(inputTracks->at(iseed));
    }
  }

  ACTS_VERBOSE("Estimated " << trackParameters.size() << " track parameters");

  m_outputTrackParameters(ctx, std::move(trackParameters));
  if (m_outputSeeds.isInitialized()) {
    m_outputSeeds(ctx, std::move(outputSeeds));
  }

  if (m_outputTracks.isInitialized()) {
    m_outputTracks(ctx, std::move(outputTracks));
  }

  return ProcessCode::SUCCESS;
}
}  // namespace ActsExamples

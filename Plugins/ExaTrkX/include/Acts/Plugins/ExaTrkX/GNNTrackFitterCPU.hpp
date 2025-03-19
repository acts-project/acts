// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/Seeding/EstimateTrackParamsFromSeed.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>

#include <boost/container/static_vector.hpp>

namespace Acts {

struct GNNParametersBuilderCPU {
  struct Config {
    std::shared_ptr<const Acts::MagneticFieldProvider> bField;
    std::shared_ptr<const Acts::TrackingGeometry> tGeometry;

    std::size_t rIdx{}, phiIdx{}, zIdx{};
    float rScale{1.f}, phiScale{1.f}, zScale{1.f};
    std::size_t nFeatures{};
    std::set<std::size_t> stripVolumes;

    EstimateTrackParamCovarianceConfig covCfg;
    ParticleHypothesis partHypot = ParticleHypothesis::pion();

    double minSpacepointDist{};
    double bFieldMin = 0.0;
    bool buildTightSeeds = true;
  };

  GNNParametersBuilderCPU(const Config &cfg,
                          std::unique_ptr<const Acts::Logger> logger)
      : m_cfg(cfg), m_logger(std::move(logger)) {
    if (!m_logger) {
      throw std::runtime_error("Missing logger!");
    }
    if (!m_cfg.bField) {
      throw std::invalid_argument("Missing bfield!");
    }
    if (!m_cfg.tGeometry) {
      throw std::invalid_argument("Missing geometry!");
    }
    if (m_cfg.nFeatures == 0) {
      throw std::invalid_argument("Cannot have 0 spacepoint features");
    }
    if (m_cfg.rIdx == m_cfg.phiIdx || m_cfg.phiIdx == m_cfg.zIdx ||
        m_cfg.zIdx == m_cfg.rIdx) {
      throw std::invalid_argument(
          "r, phi and z idx should point to different offsets!");
    }
  }

  std::optional<BoundTrackParameters> buildParameters(
      const std::vector<float> &spacepointFeatures,
      const std::vector<Acts::GeometryIdentifier> &geoIds,
      const std::vector<int> &candidate,
      Acts::MagneticFieldProvider::Cache &bCache,
      const Acts::GeometryContext &gctx) const;

 private:
  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;
  const Acts::Logger &logger() const { return *m_logger; }
};

template <typename track_container_t>
class GNNTrackFitterCPU {
 public:
  using TCBackend = typename track_container_t::TrackContainerBackend;
  using TSBackend = typename track_container_t::TrackStateContainerBackend;

  struct Config {
    std::shared_ptr<const Acts::TrackingGeometry> geometry;
    std::shared_ptr<const Acts::MagneticFieldProvider> bfield;
    GNNParametersBuilderCPU::Config paramBuilderCfg;
  };

  struct Options {
    const Acts::GeometryContext &gctx;
    const Acts::MagneticFieldContext &mctx;
    const Acts::CalibrationContext &cctx;
    KalmanFitterExtensions<TSBackend> extensions;
    const Acts::Surface *targetSurface{};
  };

  using Propagator = Acts::Propagator<Acts::SympyStepper, Acts::Navigator>;
  using Fitter = Acts::KalmanFitter<Propagator, TSBackend>;

  GNNTrackFitterCPU(const Config &cfg,
                    std::unique_ptr<const Acts::Logger> logger)
      : m_cfg(cfg), m_logger(std::move(logger)) {
    if (!m_logger) {
      throw std::invalid_argument("Missing logger!");
    }
    if (!m_cfg.geometry) {
      throw std::invalid_argument("Missing geometry!");
    }
    if (!m_cfg.bfield) {
      throw std::invalid_argument("Missing bfield!");
    }
    m_paramBuilder = std::make_unique<GNNParametersBuilderCPU>(
        cfg.paramBuilderCfg, m_logger->clone());
    Acts::Navigator::Config navCfg;
    navCfg.trackingGeometry = m_cfg.geometry;
    Acts::Navigator nav(navCfg, m_logger->clone());
    Acts::SympyStepper stepper(m_cfg.bfield);
    Propagator prop(std::move(stepper), std::move(nav), m_logger->clone());
    m_fitter = std::make_unique<Fitter>(std::move(prop), m_logger->clone());
  }

  void operator()(
      track_container_t &tracks,
      const std::vector<std::vector<int>> &candidates,
      const std::vector<float> &spacepointFeatures,
      const std::vector<Acts::GeometryIdentifier> &geoIds,
      const std::vector<boost::container::static_vector<Acts::SourceLink, 2>>
          &sourceLinks,
      Options options) const {
    auto bCache = m_cfg.bfield->makeCache(options.mctx);

    for (const auto &candidate : candidates) {
      ACTS_VERBOSE("Build parameters...");
      auto params = m_paramBuilder->buildParameters(
          spacepointFeatures, geoIds, candidate, bCache, options.gctx);

      if (!params) {
        ACTS_DEBUG("No parameters, skip candidate");
        continue;
      }

      ACTS_VERBOSE("Build options...");
      Acts::PropagatorPlainOptions popts(options.gctx, options.mctx);
      Acts::KalmanFitterOptions<TSBackend> kfOpts(
          options.gctx, options.mctx, options.cctx, options.extensions, popts);
      kfOpts.referenceSurface = options.targetSurface;
      kfOpts.referenceSurfaceStrategy =
          Acts::KalmanFitterTargetSurfaceStrategy::first;

      ACTS_VERBOSE("Collect source links...");
      std::vector<Acts::SourceLink> sls;
      for (auto i : candidate) {
        for (const auto &sl : sourceLinks.at(i)) {
          sls.push_back(sl);
        }
      }

      ACTS_VERBOSE("Start fit...");
      auto res = m_fitter->fit(sls.begin(), sls.end(), *params, kfOpts, tracks);
      if (!res.ok()) {
        ACTS_DEBUG("Track fit failed!");
        continue;
      }

      if (!res->hasReferenceSurface()) {
        ACTS_DEBUG("Fit successful, but no reference surface");
        continue;
      }
    }
  }

 private:
  Config m_cfg;
  std::unique_ptr<GNNParametersBuilderCPU> m_paramBuilder;
  std::unique_ptr<Fitter> m_fitter;

  std::unique_ptr<const Acts::Logger> m_logger;
  const Acts::Logger &logger() const { return *m_logger; }
};

}  // namespace Acts

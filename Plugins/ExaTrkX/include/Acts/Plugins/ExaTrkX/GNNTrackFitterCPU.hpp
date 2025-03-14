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
      : m_cfg(cfg), m_logger(std::move(logger)) {}

  std::optional<BoundTrackParameters> buildParameters(
      const std::vector<float> &spacepointFeatures,
      const std::vector<Acts::GeometryIdentifier> &geoIds,
      const std::vector<int> &candidate) const;

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
    KalmanFitterExtensions<TSBackend> extensions;
    GNNParametersBuilderCPU::Config paramBuilderCfg;
  };

  using Propagator = Acts::Propagator<Acts::SympyStepper, Acts::Navigator>;
  using Fitter = Acts::KalmanFitter<Propagator, TSBackend>;

  GNNTrackFitterCPU(const Config &cfg,
                    std::unique_ptr<const Acts::Logger> logger)
      : m_cfg(cfg), m_logger(std::move(logger)) {
    m_paramBuilder = std::make_unique<GNNParametersBuilderCPU>(
        cfg.paramBuilderCfg, logger->clone());
  }

  void operator()(
      track_container_t &tracks,
      const std::vector<std::vector<int>> &candidates,
      const std::vector<float> &spacepointFeatures,
      const std::vector<Acts::GeometryIdentifier> &geoIds,
      const std::vector<boost::container::static_vector<Acts::SourceLink, 2>>
          &sourceLinks,
      const Acts::GeometryContext &gctx, const Acts::MagneticFieldContext &mctx,
      const Acts::CalibrationContext &cctx) const {
    for (const auto &candidate : candidates) {
      auto params = m_paramBuilder->buildParameters(spacepointFeatures, geoIds,
                                                    candidate);

      if (!params) {
        ACTS_DEBUG("No parameters, skip candidate");
        continue;
      }

      Acts::PropagatorPlainOptions popts(gctx, mctx);
      Acts::KalmanFitterOptions<TSBackend> opts(gctx, mctx, cctx,
                                                m_cfg.extensions, popts);

      std::vector<Acts::SourceLink> sls;
      for (auto i : candidate) {
        for (const auto &sl : sourceLinks.at(i)) {
          sls.push_back(sl);
        }
      }

      auto res = m_fitter->fit(sls.begin(), sls.end(), *params, opts, tracks);
      if (!res.ok()) {
        ACTS_WARNING("Track fit failed!");
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

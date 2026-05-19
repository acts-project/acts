// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/ActsMeasToTracccAlg.hpp"

#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include <Acts/EventData/SourceLink.hpp>

#include <stdexcept>

#include <detray/geometry/identifier.hpp>

namespace ActsExamples {

ActsMeasToTracccAlg::ActsMeasToTracccAlg(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("ActsMeasToTracccAlg", std::move(logger)), m_cfg(cfg) {
  if (!m_cfg.trackingGeometry) {
    throw std::invalid_argument(
        "TracccMeasToActsAlg: trackingGeometry is null");
  }
  m_inputActsMeasurements.initialize(m_cfg.inputActsMeasurements);
  m_outputDetrayToActsMap.initialize(m_cfg.outputDetrayToActsMap);
  m_outputTracccMeasurements.initialize(m_cfg.outputTracccMeasurements);

  if (m_detrayToActsMap.empty()) {
    buildSurfaceMap(*m_cfg.trackingGeometry, m_cfg.detectorFile,
                    m_detrayToActsMap, m_actsToDetrayMap);
  }
}

void ActsExamples::ActsMeasToTracccAlg::buildSurfaceMap(
    const Acts::TrackingGeometry& trackingGeometry,
    const std::string& detectorFile,
    std::unordered_map<std::uint64_t, Acts::GeometryIdentifier>&
        detrayToActsMap,
    std::unordered_map<Acts::GeometryIdentifier, std::uint64_t>&
        actsToDetrayMap) {
  std::unordered_map<Acts::GeometryIdentifier, std::array<double, 3>>
      actsCentres;
  std::unordered_map<std::uint64_t, std::array<double, 3>> detrayCentres;

  //   Load detray detector and build Acts→detray map
  traccc::io::read_detector(m_host_det, m_mr, detectorFile, "", "");

  int nPrint = 0;
  typename traccc::default_detector::host::geometry_context detrayContext{};
  traccc::host_detector_visitor<traccc::detector_type_list>(
      m_host_det, [&]<typename detector_traits_t>(
                      const typename detector_traits_t::host& det) {
        for (const auto& surface : det.surfaces()) {
          if (!surface.is_sensitive())
            continue;
          const auto sf = detray::geometry::surface{det, surface};

          double placement_x = sf.center(detrayContext)[0];
          double placement_y = sf.center(detrayContext)[1];
          double placement_z = sf.center(detrayContext)[2];
          if (nPrint < 10) {
            ACTS_INFO("detray=" << surface.identifier().value() << " pos: "
                                << placement_x << ", " << placement_y << ", "
                                << placement_z << std::endl);
            nPrint++;
          }
          detrayCentres[surface.identifier().value()] =
              std::array{placement_x, placement_y, placement_z};
        }
      });

  ACTS_INFO("ActsMeasToTracccAlg: built detray map with "
            << detrayCentres.size() << " entries" << std::endl);

  nPrint = 0;
  Acts::GeometryContext gctx =
      Acts::GeometryContext::dangerouslyDefaultConstruct();
  trackingGeometry.visitSurfaces([&](const Acts::Surface* surface) {
    if (surface != nullptr) {
      Acts::Vector3 geo_center = surface->center(gctx);
      if (nPrint < 10) {
        ACTS_INFO("acts=" << surface->geometryId() << " pos: " << geo_center[0]
                          << ", " << geo_center[1] << ", " << geo_center[2]
                          << std::endl);
        nPrint++;
      }
      actsCentres[surface->geometryId()] =
          std::array{geo_center[0], geo_center[1], geo_center[2]};
    }
  });

  ACTS_INFO("ActsMeasToTracccAlg: built acts map with "
            << actsCentres.size() << " entries" << std::endl);

  for (const auto& [detrayID, detrayCenter] : detrayCentres) {
    for (auto& [actsID, actsCenter] : actsCentres) {
      double distance = std::sqrt(std::pow(detrayCenter[0] - actsCenter[0], 2) +
                                  std::pow(detrayCenter[1] - actsCenter[1], 2) +
                                  std::pow(detrayCenter[2] - actsCenter[2], 2));
      if (distance < 0.1) {  // threshold for matching, may need tuning
        detrayToActsMap[detrayID] = actsID;
        actsToDetrayMap[actsID] = detrayID;
        break;
      }
    }
  }

  ACTS_INFO("ActsMeasToTracccAlg: built detray to Acts map with "
            << detrayToActsMap.size() << " entries" << std::endl);
}

ProcessCode ActsMeasToTracccAlg::execute(const AlgorithmContext& ctx) const {
  const auto& actsMeasurements = m_inputActsMeasurements(ctx);
  const auto& detrayToActsMap = m_detrayToActsMap;
  // Build inverse map: Acts GeometryIdentifier → detray identifier
  m_outputDetrayToActsMap(
      ctx, std::unordered_map<std::uint64_t, Acts::GeometryIdentifier>(
               m_detrayToActsMap));
  std::unordered_map<std::uint64_t, std::uint64_t> actsToDetrayMap;
  for (const auto& [detrayId, actsId] : detrayToActsMap) {
    actsToDetrayMap[actsId.value()] = detrayId;
  }

  traccc::edm::measurement_collection::host tracccMeasurements{m_mr};

  std::size_t nConverted = 0;
  std::size_t nPix = 0, nStripShort = 0, nStripLong = 0;
  for (std::size_t i = 0; i < actsMeasurements.size(); ++i) {
    const auto meas = actsMeasurements.getMeasurement(i);
    const auto geoId = meas.geometryId();
    const auto dims = static_cast<unsigned int>(meas.size());

    if (dims == 0u)
      continue;

    auto it = actsToDetrayMap.find(geoId.value());
    if (it == actsToDetrayMap.end()) {
      ACTS_WARNING("ActsMeasToTracccAlg: no detray id for Acts geo id "
                   << geoId);
      continue;
    }

    tracccMeasurements.push_back({});
    auto tm = tracccMeasurements.at(tracccMeasurements.size() - 1);

    const auto vol = geoId.volume();
    if (std::find(m_cfg.pixelVolumes.begin(), m_cfg.pixelVolumes.end(), vol) !=
        m_cfg.pixelVolumes.end()) {
      nPix++;
      // Pixel (dims=2 or 3 with timing)
      tm.dimensions() = 2u;
      tm.local_position()[0] = meas.parameters()[Acts::eBoundLoc0];
      tm.local_position()[1] = meas.parameters()[Acts::eBoundLoc1];
      tm.local_variance()[0] = meas.covariance()(0, 0);
      tm.local_variance()[1] = meas.covariance()(1, 1);
      tm.subspace()[0] = 0u;
      tm.subspace()[1] = 1u;

    } else {
      // 1D measurement - determine if long strip (1D in traccc) or short strip
      // (2D in traccc)

      if (dims == 2u) {
        // Short strip barrel: 2D in traccc with large cross-strip variance
        tm.dimensions() = 2u;
        tm.local_position()[0] = meas.parameters()[Acts::eBoundLoc0];
        tm.local_position()[1] = meas.parameters()[Acts::eBoundLoc1];
        tm.local_variance()[0] = meas.covariance()(0, 0);
        tm.local_variance()[1] = meas.covariance()(1, 1);
        tm.subspace()[0] = 0u;
        tm.subspace()[1] = 1u;
        nStripShort++;

      } else {
        // Everything else (long strips barrel/endcap, short strip endcap): 1D
        tm.dimensions() = 1u;
        tm.local_position()[0] = meas.parameters()[Acts::eBoundLoc0];
        tm.local_position()[1] = 0.f;
        tm.local_variance()[0] = meas.covariance()(0, 0);
        tm.local_variance()[1] = 0.f;
        tm.subspace()[0] = 0u;
        tm.subspace()[1] = 1u;
        nStripLong++;
      }
    }

    tm.surface_link() = detray::geometry::identifier(it->second);
    tm.time() = 0.f;
    tm.diameter() = 0.f;
    tm.identifier() = static_cast<unsigned int>(i);  // store Acts index
    tm.cluster_index() = 0u;

    ++nConverted;
  }

  ACTS_INFO("ActsMeasToTracccAlg: converted " << nConverted << " / "
                                              << actsMeasurements.size()
                                              << " measurements");
  ACTS_INFO("  of which " << nPix << " pixel, " << nStripShort
                          << " short strip, " << nStripLong << " long strip");

  ACTS_INFO("ActsMeasToTracccAlg: built " << tracccMeasurements.size()
                                          << " measurements");

  m_outputTracccMeasurements(ctx, std::move(tracccMeasurements));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

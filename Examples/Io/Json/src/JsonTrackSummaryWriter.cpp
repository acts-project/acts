// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Json/JsonTrackSummaryWriter.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/AnyTrackProxy.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/GsfOptions.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/detail/periodic.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <numbers>
#include <optional>
#include <stdexcept>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

namespace ActsExamples {

JsonTrackSummaryWriter::JsonTrackSummaryWriter(
    const JsonTrackSummaryWriter::Config& config, Acts::Logging::Level level)
    : WriterT(config.inputTracks, "JsonTrackSummaryWriter", level),
      m_cfg(config) {
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }

  m_inputParticles.maybeInitialize(m_cfg.inputParticles);
  m_inputTrackParticleMatching.maybeInitialize(
      m_cfg.inputTrackParticleMatching);
}

JsonTrackSummaryWriter::~JsonTrackSummaryWriter() = default;

ProcessCode JsonTrackSummaryWriter::finalize() {
  nlohmann::json root = m_events;

  std::ofstream ofstream(m_cfg.filePath);
  if (!ofstream.is_open()) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }
  ofstream << root.dump(2);
  ACTS_INFO("Wrote track summary to '" << m_cfg.filePath << "'");

  return ProcessCode::SUCCESS;
}

ProcessCode JsonTrackSummaryWriter::writeT(const AlgorithmContext& ctx,
                                           const ConstTrackContainer& tracks) {
  const static SimParticleContainer emptyParticles;
  const static TrackParticleMatching emptyTrackParticleMatching;

  const auto& particles =
      m_inputParticles.isInitialized() ? m_inputParticles(ctx) : emptyParticles;
  const auto& trackParticleMatching =
      m_inputTrackParticleMatching.isInitialized()
          ? m_inputTrackParticleMatching(ctx)
          : emptyTrackParticleMatching;

  std::vector<ParticleHitCount> particleHitCounts;

  nlohmann::json tracksJson = nlohmann::json::array();

  for (const auto& track : tracks) {
    nlohmann::json t;

    t["track_nr"] = track.index();
    t["nStates"] = track.nTrackStates();
    t["nMeasurements"] = track.nMeasurements();
    t["nOutliers"] = track.nOutliers();
    t["nHoles"] = track.nHoles();
    t["nSharedHits"] = track.nSharedHits();
    t["chi2Sum"] = track.chi2();
    t["NDF"] = track.nDoF();

    // Per-state chi2 and geometry IDs
    {
      std::vector<double> mChi2, oChi2;
      std::vector<std::uint32_t> mVol, mLay, oVol, oLay;
      for (const auto& state : track.trackStatesReversed()) {
        const auto& geoID = state.referenceSurface().geometryId();
        if (state.typeFlags().isOutlier()) {
          oChi2.push_back(state.chi2());
          oVol.push_back(geoID.volume());
          oLay.push_back(geoID.layer());
        } else if (state.typeFlags().isMeasurement()) {
          mChi2.push_back(state.chi2());
          mVol.push_back(geoID.volume());
          mLay.push_back(geoID.layer());
        }
      }
      t["measurementChi2"] = std::move(mChi2);
      t["outlierChi2"] = std::move(oChi2);
      t["measurementVolume"] = std::move(mVol);
      t["measurementLayer"] = std::move(mLay);
      t["outlierVolume"] = std::move(oVol);
      t["outlierLayer"] = std::move(oLay);
    }

    // Truth particle matching
    SimBarcode majorityParticleId{};
    TrackMatchClassification trackClassification =
        TrackMatchClassification::Unknown;
    unsigned int nMajorityHits = std::numeric_limits<unsigned int>::max();
    int t_charge = std::numeric_limits<int>::max();
    float t_time = NaNfloat;
    float t_vx = NaNfloat, t_vy = NaNfloat, t_vz = NaNfloat;
    float t_px = NaNfloat, t_py = NaNfloat, t_pz = NaNfloat;
    float t_theta_v = NaNfloat, t_phi_v = NaNfloat, t_eta_v = NaNfloat;
    float t_p = NaNfloat, t_pT = NaNfloat;
    float t_d0 = NaNfloat, t_z0 = NaNfloat, t_prodR = NaNfloat;
    float t_qop = NaNfloat;

    const Acts::Surface* pSurface =
        track.hasReferenceSurface() ? &track.referenceSurface() : nullptr;

    bool foundMajorityParticle = false;
    auto match = trackParticleMatching.find(track.index());
    if (match != trackParticleMatching.end() &&
        match->second.particle.has_value()) {
      majorityParticleId = match->second.particle.value();
      trackClassification = match->second.classification;
      nMajorityHits = match->second.contributingParticles.front().hitCount;

      auto ip = particles.find(majorityParticleId);
      if (ip != particles.end()) {
        foundMajorityParticle = true;
        const auto& p = *ip;
        t_p = p.absoluteMomentum();
        t_charge = static_cast<int>(p.charge());
        t_time = p.time();
        t_vx = p.position().x();
        t_vy = p.position().y();
        t_vz = p.position().z();
        t_px = t_p * p.direction().x();
        t_py = t_p * p.direction().y();
        t_pz = t_p * p.direction().z();
        t_theta_v = theta(p.direction());
        t_phi_v = phi(p.direction());
        t_eta_v = eta(p.direction());
        t_pT = t_p * perp(p.direction());
        t_qop = p.qOverP();
        t_prodR = std::sqrt(t_vx * t_vx + t_vy * t_vy);

        if (pSurface != nullptr) {
          Acts::Intersection3D intersection =
              pSurface
                  ->intersect(ctx.geoContext, p.position(), p.direction(),
                              Acts::BoundaryTolerance::Infinite())
                  .closest();
          auto pos = intersection.position();
          auto lpResult =
              pSurface->globalToLocal(ctx.geoContext, pos, p.direction());
          if (lpResult.ok()) {
            t_d0 = lpResult.value()[Acts::BoundIndices::eBoundLoc0];
            t_z0 = lpResult.value()[Acts::BoundIndices::eBoundLoc1];
          } else {
            ACTS_ERROR("Global to local transformation did not succeed.");
          }
        }
      }
    }

    t["nMajorityHits"] = nMajorityHits;
    t["majorityParticleId_vertex_primary"] = majorityParticleId.vertexPrimary();
    t["majorityParticleId_vertex_secondary"] =
        majorityParticleId.vertexSecondary();
    t["majorityParticleId_particle"] = majorityParticleId.particle();
    t["majorityParticleId_generation"] = majorityParticleId.generation();
    t["majorityParticleId_sub_particle"] = majorityParticleId.subParticle();
    t["trackClassification"] = static_cast<int>(trackClassification);
    t["t_charge"] = t_charge;
    t["t_time"] = t_time;
    t["t_vx"] = t_vx;
    t["t_vy"] = t_vy;
    t["t_vz"] = t_vz;
    t["t_px"] = t_px;
    t["t_py"] = t_py;
    t["t_pz"] = t_pz;
    t["t_theta"] = t_theta_v;
    t["t_phi"] = t_phi_v;
    t["t_eta"] = t_eta_v;
    t["t_p"] = t_p;
    t["t_pT"] = t_pT;
    t["t_d0"] = t_d0;
    t["t_z0"] = t_z0;
    t["t_prodR"] = t_prodR;

    // Fitted track parameters
    std::array<float, Acts::eBoundSize> param = {NaNfloat, NaNfloat, NaNfloat,
                                                 NaNfloat, NaNfloat, NaNfloat};
    std::array<float, Acts::eBoundSize> error = {NaNfloat, NaNfloat, NaNfloat,
                                                 NaNfloat, NaNfloat, NaNfloat};

    bool hasFittedParams = track.hasReferenceSurface();
    if (hasFittedParams) {
      const auto& parameter = track.parameters();
      for (unsigned int i = 0; i < Acts::eBoundSize; ++i) {
        param[i] = parameter[i];
      }
      for (unsigned int i = 0; i < Acts::eBoundSize; ++i) {
        double variance = track.covariance()(i, i);
        error[i] = variance >= 0 ? std::sqrt(variance) : NaNfloat;
      }
    }

    std::array<float, Acts::eBoundSize> res = {NaNfloat, NaNfloat, NaNfloat,
                                               NaNfloat, NaNfloat, NaNfloat};
    std::array<float, Acts::eBoundSize> pull = {NaNfloat, NaNfloat, NaNfloat,
                                                NaNfloat, NaNfloat, NaNfloat};
    if (foundMajorityParticle && hasFittedParams) {
      res = {param[Acts::eBoundLoc0] - t_d0,
             param[Acts::eBoundLoc1] - t_z0,
             Acts::detail::difference_periodic(
                 param[Acts::eBoundPhi], t_phi_v,
                 static_cast<float>(2 * std::numbers::pi)),
             param[Acts::eBoundTheta] - t_theta_v,
             param[Acts::eBoundQOverP] - t_qop,
             param[Acts::eBoundTime] - t_time};
      for (unsigned int i = 0; i < Acts::eBoundSize; ++i) {
        pull[i] = res[i] / error[i];
      }
    }

    t["hasFittedParams"] = hasFittedParams;
    t["eLOC0_fit"] = param[Acts::eBoundLoc0];
    t["eLOC1_fit"] = param[Acts::eBoundLoc1];
    t["ePHI_fit"] = param[Acts::eBoundPhi];
    t["eTHETA_fit"] = param[Acts::eBoundTheta];
    t["eQOP_fit"] = param[Acts::eBoundQOverP];
    t["eT_fit"] = param[Acts::eBoundTime];

    t["err_eLOC0_fit"] = error[Acts::eBoundLoc0];
    t["err_eLOC1_fit"] = error[Acts::eBoundLoc1];
    t["err_ePHI_fit"] = error[Acts::eBoundPhi];
    t["err_eTHETA_fit"] = error[Acts::eBoundTheta];
    t["err_eQOP_fit"] = error[Acts::eBoundQOverP];
    t["err_eT_fit"] = error[Acts::eBoundTime];

    t["res_eLOC0_fit"] = res[Acts::eBoundLoc0];
    t["res_eLOC1_fit"] = res[Acts::eBoundLoc1];
    t["res_ePHI_fit"] = res[Acts::eBoundPhi];
    t["res_eTHETA_fit"] = res[Acts::eBoundTheta];
    t["res_eQOP_fit"] = res[Acts::eBoundQOverP];
    t["res_eT_fit"] = res[Acts::eBoundTime];

    t["pull_eLOC0_fit"] = pull[Acts::eBoundLoc0];
    t["pull_eLOC1_fit"] = pull[Acts::eBoundLoc1];
    t["pull_ePHI_fit"] = pull[Acts::eBoundPhi];
    t["pull_eTHETA_fit"] = pull[Acts::eBoundTheta];
    t["pull_eQOP_fit"] = pull[Acts::eBoundQOverP];
    t["pull_eT_fit"] = pull[Acts::eBoundTime];

    // Optional covariance matrix
    if (m_cfg.writeCovMat) {
      nlohmann::json cov;
      for (int i = 0; i < static_cast<int>(Acts::eBoundSize); ++i) {
        for (int j = 0; j < static_cast<int>(Acts::eBoundSize); ++j) {
          std::string key =
              "cov_" + std::to_string(i) + "_" + std::to_string(j);
          cov[key] = track.covariance()(i, j);
        }
      }
      t["covariance"] = std::move(cov);
    }

    // Optional GSF material statistics
    if (m_cfg.writeGsfSpecific) {
      using namespace Acts::GsfConstants;
      float maxMat = NaNfloat, sumMat = NaNfloat;
      if (tracks.hasColumn(Acts::hashString(kFwdMaxMaterialXOverX0))) {
        maxMat = static_cast<float>(
            track.template component<double>(kFwdMaxMaterialXOverX0));
      }
      if (tracks.hasColumn(Acts::hashString(kFwdSumMaterialXOverX0))) {
        sumMat = static_cast<float>(
            track.template component<double>(kFwdSumMaterialXOverX0));
      }
      t["gsf_max_material_fwd"] = maxMat;
      t["gsf_sum_material_fwd"] = sumMat;
    }

    // Optional GX2F update count
    if (m_cfg.writeGx2fSpecific) {
      int nUpdate = -1;
      if (tracks.hasColumn(Acts::hashString("Gx2fnUpdateColumn"))) {
        nUpdate = static_cast<int>(
            track.template component<std::uint32_t,
                                     Acts::hashString("Gx2fnUpdateColumn")>());
      }
      t["nUpdatesGx2f"] = nUpdate;
    }

    tracksJson.push_back(std::move(t));
  }

  nlohmann::json eventJson;
  eventJson["event_nr"] = ctx.eventNumber;
  eventJson["tracks"] = std::move(tracksJson);

  std::lock_guard<std::mutex> lock(m_writeMutex);
  m_events.push_back(std::move(eventJson));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

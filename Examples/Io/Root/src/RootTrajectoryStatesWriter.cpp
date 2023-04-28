// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootTrajectoryStatesWriter.hpp"

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/detail/periodic.hpp"
#include "ActsExamples/EventData/AverageSimHits.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Utilities/Range.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <ios>
#include <limits>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

ActsExamples::RootTrajectoryStatesWriter::RootTrajectoryStatesWriter(
    const ActsExamples::RootTrajectoryStatesWriter::Config& config,
    Acts::Logging::Level level)
    : WriterT(config.inputTrajectories, "RootTrajectoryStatesWriter", level),
      m_cfg(config) {
  // trajectories collection name is already checked by base ctor
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing particles input collection");
  }
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing hit-particles map input collection");
  }
  if (m_cfg.inputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-simulated-hits map input collection");
  }
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputMeasurementParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);
  m_inputMeasurementSimHitsMap.initialize(m_cfg.inputMeasurementSimHitsMap);

  // Setup ROOT I/O
  auto path = m_cfg.filePath;
  m_outputFile = TFile::Open(path.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + path);
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  } else {
    // I/O parameters
    m_outputTree->Branch("event_nr", &m_eventNr);
    m_outputTree->Branch("multiTraj_nr", &m_multiTrajNr);
    m_outputTree->Branch("subTraj_nr", &m_subTrajNr);

    m_outputTree->Branch("t_x", &m_t_x);
    m_outputTree->Branch("t_y", &m_t_y);
    m_outputTree->Branch("t_z", &m_t_z);
    m_outputTree->Branch("t_r", &m_t_r);
    m_outputTree->Branch("t_dx", &m_t_dx);
    m_outputTree->Branch("t_dy", &m_t_dy);
    m_outputTree->Branch("t_dz", &m_t_dz);
    m_outputTree->Branch("t_eLOC0", &m_t_eLOC0);
    m_outputTree->Branch("t_eLOC1", &m_t_eLOC1);
    m_outputTree->Branch("t_ePHI", &m_t_ePHI);
    m_outputTree->Branch("t_eTHETA", &m_t_eTHETA);
    m_outputTree->Branch("t_eQOP", &m_t_eQOP);
    m_outputTree->Branch("t_eT", &m_t_eT);

    m_outputTree->Branch("nStates", &m_nStates);
    m_outputTree->Branch("nMeasurements", &m_nMeasurements);
    m_outputTree->Branch("volume_id", &m_volumeID);
    m_outputTree->Branch("layer_id", &m_layerID);
    m_outputTree->Branch("module_id", &m_moduleID);
    m_outputTree->Branch("pathLength", &m_pathLength);
    m_outputTree->Branch("l_x_hit", &m_lx_hit);
    m_outputTree->Branch("l_y_hit", &m_ly_hit);
    m_outputTree->Branch("g_x_hit", &m_x_hit);
    m_outputTree->Branch("g_y_hit", &m_y_hit);
    m_outputTree->Branch("g_z_hit", &m_z_hit);
    m_outputTree->Branch("res_x_hit", &m_res_x_hit);
    m_outputTree->Branch("res_y_hit", &m_res_y_hit);
    m_outputTree->Branch("err_x_hit", &m_err_x_hit);
    m_outputTree->Branch("err_y_hit", &m_err_y_hit);
    m_outputTree->Branch("pull_x_hit", &m_pull_x_hit);
    m_outputTree->Branch("pull_y_hit", &m_pull_y_hit);
    m_outputTree->Branch("dim_hit", &m_dim_hit);

    m_outputTree->Branch("nPredicted", &m_nParams[0]);
    m_outputTree->Branch("predicted", &m_hasParams[0]);
    m_outputTree->Branch("eLOC0_prt", &m_eLOC0[0]);
    m_outputTree->Branch("eLOC1_prt", &m_eLOC1[0]);
    m_outputTree->Branch("ePHI_prt", &m_ePHI[0]);
    m_outputTree->Branch("eTHETA_prt", &m_eTHETA[0]);
    m_outputTree->Branch("eQOP_prt", &m_eQOP[0]);
    m_outputTree->Branch("eT_prt", &m_eT[0]);
    m_outputTree->Branch("res_eLOC0_prt", &m_res_eLOC0[0]);
    m_outputTree->Branch("res_eLOC1_prt", &m_res_eLOC1[0]);
    m_outputTree->Branch("res_ePHI_prt", &m_res_ePHI[0]);
    m_outputTree->Branch("res_eTHETA_prt", &m_res_eTHETA[0]);
    m_outputTree->Branch("res_eQOP_prt", &m_res_eQOP[0]);
    m_outputTree->Branch("res_eT_prt", &m_res_eT[0]);
    m_outputTree->Branch("err_eLOC0_prt", &m_err_eLOC0[0]);
    m_outputTree->Branch("err_eLOC1_prt", &m_err_eLOC1[0]);
    m_outputTree->Branch("err_ePHI_prt", &m_err_ePHI[0]);
    m_outputTree->Branch("err_eTHETA_prt", &m_err_eTHETA[0]);
    m_outputTree->Branch("err_eQOP_prt", &m_err_eQOP[0]);
    m_outputTree->Branch("err_eT_prt", &m_err_eT[0]);
    m_outputTree->Branch("pull_eLOC0_prt", &m_pull_eLOC0[0]);
    m_outputTree->Branch("pull_eLOC1_prt", &m_pull_eLOC1[0]);
    m_outputTree->Branch("pull_ePHI_prt", &m_pull_ePHI[0]);
    m_outputTree->Branch("pull_eTHETA_prt", &m_pull_eTHETA[0]);
    m_outputTree->Branch("pull_eQOP_prt", &m_pull_eQOP[0]);
    m_outputTree->Branch("pull_eT_prt", &m_pull_eT[0]);
    m_outputTree->Branch("g_x_prt", &m_x[0]);
    m_outputTree->Branch("g_y_prt", &m_y[0]);
    m_outputTree->Branch("g_z_prt", &m_z[0]);
    m_outputTree->Branch("px_prt", &m_px[0]);
    m_outputTree->Branch("py_prt", &m_py[0]);
    m_outputTree->Branch("pz_prt", &m_pz[0]);
    m_outputTree->Branch("eta_prt", &m_eta[0]);
    m_outputTree->Branch("pT_prt", &m_pT[0]);

    m_outputTree->Branch("nFiltered", &m_nParams[1]);
    m_outputTree->Branch("filtered", &m_hasParams[1]);
    m_outputTree->Branch("eLOC0_flt", &m_eLOC0[1]);
    m_outputTree->Branch("eLOC1_flt", &m_eLOC1[1]);
    m_outputTree->Branch("ePHI_flt", &m_ePHI[1]);
    m_outputTree->Branch("eTHETA_flt", &m_eTHETA[1]);
    m_outputTree->Branch("eQOP_flt", &m_eQOP[1]);
    m_outputTree->Branch("eT_flt", &m_eT[1]);
    m_outputTree->Branch("res_eLOC0_flt", &m_res_eLOC0[1]);
    m_outputTree->Branch("res_eLOC1_flt", &m_res_eLOC1[1]);
    m_outputTree->Branch("res_ePHI_flt", &m_res_ePHI[1]);
    m_outputTree->Branch("res_eTHETA_flt", &m_res_eTHETA[1]);
    m_outputTree->Branch("res_eQOP_flt", &m_res_eQOP[1]);
    m_outputTree->Branch("res_eT_flt", &m_res_eT[1]);
    m_outputTree->Branch("err_eLOC0_flt", &m_err_eLOC0[1]);
    m_outputTree->Branch("err_eLOC1_flt", &m_err_eLOC1[1]);
    m_outputTree->Branch("err_ePHI_flt", &m_err_ePHI[1]);
    m_outputTree->Branch("err_eTHETA_flt", &m_err_eTHETA[1]);
    m_outputTree->Branch("err_eQOP_flt", &m_err_eQOP[1]);
    m_outputTree->Branch("err_eT_flt", &m_err_eT[1]);
    m_outputTree->Branch("pull_eLOC0_flt", &m_pull_eLOC0[1]);
    m_outputTree->Branch("pull_eLOC1_flt", &m_pull_eLOC1[1]);
    m_outputTree->Branch("pull_ePHI_flt", &m_pull_ePHI[1]);
    m_outputTree->Branch("pull_eTHETA_flt", &m_pull_eTHETA[1]);
    m_outputTree->Branch("pull_eQOP_flt", &m_pull_eQOP[1]);
    m_outputTree->Branch("pull_eT_flt", &m_pull_eT[1]);
    m_outputTree->Branch("g_x_flt", &m_x[1]);
    m_outputTree->Branch("g_y_flt", &m_y[1]);
    m_outputTree->Branch("g_z_flt", &m_z[1]);
    m_outputTree->Branch("px_flt", &m_px[1]);
    m_outputTree->Branch("py_flt", &m_py[1]);
    m_outputTree->Branch("pz_flt", &m_pz[1]);
    m_outputTree->Branch("eta_flt", &m_eta[1]);
    m_outputTree->Branch("pT_flt", &m_pT[1]);

    m_outputTree->Branch("nSmoothed", &m_nParams[2]);
    m_outputTree->Branch("smoothed", &m_hasParams[2]);
    m_outputTree->Branch("eLOC0_smt", &m_eLOC0[2]);
    m_outputTree->Branch("eLOC1_smt", &m_eLOC1[2]);
    m_outputTree->Branch("ePHI_smt", &m_ePHI[2]);
    m_outputTree->Branch("eTHETA_smt", &m_eTHETA[2]);
    m_outputTree->Branch("eQOP_smt", &m_eQOP[2]);
    m_outputTree->Branch("eT_smt", &m_eT[2]);
    m_outputTree->Branch("res_eLOC0_smt", &m_res_eLOC0[2]);
    m_outputTree->Branch("res_eLOC1_smt", &m_res_eLOC1[2]);
    m_outputTree->Branch("res_ePHI_smt", &m_res_ePHI[2]);
    m_outputTree->Branch("res_eTHETA_smt", &m_res_eTHETA[2]);
    m_outputTree->Branch("res_eQOP_smt", &m_res_eQOP[2]);
    m_outputTree->Branch("res_eT_smt", &m_res_eT[2]);
    m_outputTree->Branch("err_eLOC0_smt", &m_err_eLOC0[2]);
    m_outputTree->Branch("err_eLOC1_smt", &m_err_eLOC1[2]);
    m_outputTree->Branch("err_ePHI_smt", &m_err_ePHI[2]);
    m_outputTree->Branch("err_eTHETA_smt", &m_err_eTHETA[2]);
    m_outputTree->Branch("err_eQOP_smt", &m_err_eQOP[2]);
    m_outputTree->Branch("err_eT_smt", &m_err_eT[2]);
    m_outputTree->Branch("pull_eLOC0_smt", &m_pull_eLOC0[2]);
    m_outputTree->Branch("pull_eLOC1_smt", &m_pull_eLOC1[2]);
    m_outputTree->Branch("pull_ePHI_smt", &m_pull_ePHI[2]);
    m_outputTree->Branch("pull_eTHETA_smt", &m_pull_eTHETA[2]);
    m_outputTree->Branch("pull_eQOP_smt", &m_pull_eQOP[2]);
    m_outputTree->Branch("pull_eT_smt", &m_pull_eT[2]);
    m_outputTree->Branch("g_x_smt", &m_x[2]);
    m_outputTree->Branch("g_y_smt", &m_y[2]);
    m_outputTree->Branch("g_z_smt", &m_z[2]);
    m_outputTree->Branch("px_smt", &m_px[2]);
    m_outputTree->Branch("py_smt", &m_py[2]);
    m_outputTree->Branch("pz_smt", &m_pz[2]);
    m_outputTree->Branch("eta_smt", &m_eta[2]);
    m_outputTree->Branch("pT_smt", &m_pT[2]);

    m_outputTree->Branch("chi2", &m_chi2);
  }
}

ActsExamples::RootTrajectoryStatesWriter::~RootTrajectoryStatesWriter() {
  m_outputFile->Close();
}

ActsExamples::ProcessCode ActsExamples::RootTrajectoryStatesWriter::finalize() {
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();

  ACTS_INFO("Wrote states of trajectories to tree '"
            << m_cfg.treeName << "' in '" << m_cfg.treeName << "'");

  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::RootTrajectoryStatesWriter::writeT(
    const AlgorithmContext& ctx, const TrajectoriesContainer& trajectories) {
  auto& gctx = ctx.geoContext;
  // Read additional input collections
  const auto& particles = m_inputParticles(ctx);
  const auto& simHits = m_inputSimHits(ctx);
  const auto& hitParticlesMap = m_inputMeasurementParticlesMap(ctx);
  const auto& hitSimHitsMap = m_inputMeasurementSimHitsMap(ctx);

  // For each particle within a track, how many hits did it contribute
  std::vector<ParticleHitCount> particleHitCounts;

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Get the event number
  m_eventNr = ctx.eventNumber;

  // Loop over the trajectories
  for (size_t itraj = 0; itraj < trajectories.size(); ++itraj) {
    const auto& traj = trajectories[itraj];

    // The trajectory index
    m_multiTrajNr = itraj;

    // The trajectory entry indices
    const auto& trackTips = traj.tips();

    // Dont write empty MultiTrajectory
    if (trackTips.empty()) {
      continue;
    }

    // The MultiTrajectory
    const auto& mj = traj.multiTrajectory();

    // Loop over the entry indices for the subtrajectories
    for (unsigned int isubtraj = 0; isubtraj < trackTips.size(); ++isubtraj) {
      // The subtrajectory index
      m_subTrajNr = isubtraj;
      // The entry index for this subtrajectory
      const auto& trackTip = trackTips[isubtraj];
      // Collect the trajectory summary info
      auto trajState =
          Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
      m_nMeasurements = trajState.nMeasurements;
      m_nStates = trajState.nStates;

      // Get the majority truth particle to this track
      int truthQ = 1.;
      identifyContributingParticles(hitParticlesMap, traj, trackTip,
                                    particleHitCounts);
      if (not particleHitCounts.empty()) {
        // Get the barcode of the majority truth particle
        auto barcode = particleHitCounts.front().particleId;
        // Find the truth particle via the barcode
        auto ip = particles.find(barcode);
        if (ip != particles.end()) {
          const auto& particle = *ip;
          ACTS_VERBOSE("Find the truth particle with barcode "
                       << barcode << "=" << barcode.value());
          // Get the truth particle charge
          truthQ = static_cast<int>(particle.charge());
        } else {
          ACTS_DEBUG("Truth particle with barcode "
                     << barcode << "=" << barcode.value() << " not found!");
        }
      }

      // Get the trackStates on the trajectory
      m_nParams = {0, 0, 0};
      mj.visitBackwards(trackTip, [&](const auto& state) {
        // we only fill the track states with non-outlier measurement
        auto typeFlags = state.typeFlags();
        if (not typeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
          return true;
        }

        const auto& surface = state.referenceSurface();

        // get the truth hits corresponding to this trackState
        // Use average truth in the case of multiple contributing sim hits
        auto sl =
            state.getUncalibratedSourceLink().template get<IndexSourceLink>();
        const auto hitIdx = sl.index();
        auto indices = makeRange(hitSimHitsMap.equal_range(hitIdx));
        auto [truthLocal, truthPos4, truthUnitDir] =
            averageSimHits(ctx.geoContext, surface, simHits, indices);
        // momemtum averaging makes even less sense than averaging position and
        // direction. use the first momentum or set q/p to zero
        float truthQOP = 0.0f;
        if (not indices.empty()) {
          // we assume that the indices are within valid ranges so we do not
          // need to check their validity again.
          const auto simHitIdx0 = indices.begin()->second;
          const auto& simHit0 = *simHits.nth(simHitIdx0);
          const auto p =
              simHit0.momentum4Before().template segment<3>(Acts::eMom0).norm();
          truthQOP = truthQ / p;
        }

        // fill the truth hit info
        m_t_x.push_back(truthPos4[Acts::ePos0]);
        m_t_y.push_back(truthPos4[Acts::ePos1]);
        m_t_z.push_back(truthPos4[Acts::ePos2]);
        m_t_r.push_back(perp(truthPos4.template segment<3>(Acts::ePos0)));
        m_t_dx.push_back(truthUnitDir[Acts::eMom0]);
        m_t_dy.push_back(truthUnitDir[Acts::eMom1]);
        m_t_dz.push_back(truthUnitDir[Acts::eMom2]);

        // get the truth track parameter at this track State
        float truthLOC0 = truthLocal[Acts::ePos0];
        float truthLOC1 = truthLocal[Acts::ePos1];
        float truthTIME = truthPos4[Acts::eTime];
        float truthPHI = phi(truthUnitDir);
        float truthTHETA = theta(truthUnitDir);

        // fill the truth track parameter at this track State
        m_t_eLOC0.push_back(truthLOC0);
        m_t_eLOC1.push_back(truthLOC1);
        m_t_ePHI.push_back(truthPHI);
        m_t_eTHETA.push_back(truthTHETA);
        m_t_eQOP.push_back(truthQOP);
        m_t_eT.push_back(truthTIME);

        // get the geometry ID
        auto geoID = surface.geometryId();
        m_volumeID.push_back(geoID.volume());
        m_layerID.push_back(geoID.layer());
        m_moduleID.push_back(geoID.sensitive());

        // get the path length
        m_pathLength.push_back(state.pathLength());

        // expand the local measurements into the full bound space
        Acts::BoundVector meas = state.effectiveProjector().transpose() *
                                 state.effectiveCalibrated();
        // extract local and global position
        Acts::Vector2 local(meas[Acts::eBoundLoc0], meas[Acts::eBoundLoc1]);
        Acts::Vector3 mom(1, 1, 1);
        Acts::Vector3 global =
            surface.localToGlobal(ctx.geoContext, local, mom);

        // fill the measurement info
        m_lx_hit.push_back(local[Acts::ePos0]);
        m_ly_hit.push_back(local[Acts::ePos1]);
        m_x_hit.push_back(global[Acts::ePos0]);
        m_y_hit.push_back(global[Acts::ePos1]);
        m_z_hit.push_back(global[Acts::ePos2]);

        // status of the fitted track parameters
        std::array<bool, 3> hasParams = {false, false, false};
        // optional fitted track parameters
        std::optional<std::pair<Acts::BoundVector, Acts::BoundMatrix>>
            trackParamsOpt = std::nullopt;
        // lambda to get the fitted track parameters
        auto getTrackParams = [&](unsigned int ipar) {
          if (ipar == 0 && state.hasPredicted()) {
            hasParams[0] = true;
            m_nParams[0]++;
            trackParamsOpt =
                std::make_pair(state.predicted(), state.predictedCovariance());
          } else if (ipar == 1 && state.hasFiltered()) {
            hasParams[1] = true;
            m_nParams[1]++;
            trackParamsOpt =
                std::make_pair(state.filtered(), state.filteredCovariance());
          } else if (ipar == 2 && state.hasSmoothed()) {
            hasParams[2] = true;
            m_nParams[2]++;
            trackParamsOpt =
                std::make_pair(state.smoothed(), state.smoothedCovariance());
          }
        };

        // fill the fitted track parameters
        for (unsigned int ipar = 0; ipar < 3; ++ipar) {
          // get the fitted track parameters
          getTrackParams(ipar);
          if (trackParamsOpt) {
            const auto& [parameters, covariance] = *trackParamsOpt;
            if (ipar == 0) {
              //
              // local hit residual info
              auto H = state.effectiveProjector();
              auto hitCov = state.effectiveCalibratedCovariance();
              auto resCov = hitCov + H * covariance * H.transpose();
              auto res = state.effectiveCalibrated() - H * parameters;

              m_res_x_hit.push_back(res[Acts::eBoundLoc0]);
              m_err_x_hit.push_back(
                  sqrt(hitCov(Acts::eBoundLoc0, Acts::eBoundLoc0)));
              m_pull_x_hit.push_back(
                  res[Acts::eBoundLoc0] /
                  sqrt(resCov(Acts::eBoundLoc0, Acts::eBoundLoc0)));

              if (state.calibratedSize() >= 2) {
                m_pull_y_hit.push_back(
                    res[Acts::eBoundLoc1] /
                    sqrt(resCov(Acts::eBoundLoc1, Acts::eBoundLoc1)));
                m_res_y_hit.push_back(res[Acts::eBoundLoc1]);
                m_err_y_hit.push_back(
                    sqrt(hitCov(Acts::eBoundLoc1, Acts::eBoundLoc1)));
              } else {
                float nan = std::numeric_limits<float>::quiet_NaN();
                m_pull_y_hit.push_back(nan);
                m_res_y_hit.push_back(nan);
                m_err_y_hit.push_back(nan);
              }

              m_dim_hit.push_back(state.calibratedSize());
            }

            // track parameters
            m_eLOC0[ipar].push_back(parameters[Acts::eBoundLoc0]);
            m_eLOC1[ipar].push_back(parameters[Acts::eBoundLoc1]);
            m_ePHI[ipar].push_back(parameters[Acts::eBoundPhi]);
            m_eTHETA[ipar].push_back(parameters[Acts::eBoundTheta]);
            m_eQOP[ipar].push_back(parameters[Acts::eBoundQOverP]);
            m_eT[ipar].push_back(parameters[Acts::eBoundTime]);

            // track parameters residual
            m_res_eLOC0[ipar].push_back(parameters[Acts::eBoundLoc0] -
                                        truthLOC0);
            m_res_eLOC1[ipar].push_back(parameters[Acts::eBoundLoc1] -
                                        truthLOC1);
            float resPhi = Acts::detail::difference_periodic<float>(
                parameters[Acts::eBoundPhi], truthPHI,
                static_cast<float>(2 * M_PI));
            m_res_ePHI[ipar].push_back(resPhi);
            m_res_eTHETA[ipar].push_back(parameters[Acts::eBoundTheta] -
                                         truthTHETA);
            m_res_eQOP[ipar].push_back(parameters[Acts::eBoundQOverP] -
                                       truthQOP);
            m_res_eT[ipar].push_back(parameters[Acts::eBoundTime] - truthTIME);

            // track parameters error
            m_err_eLOC0[ipar].push_back(
                sqrt(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)));
            m_err_eLOC1[ipar].push_back(
                sqrt(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)));
            m_err_ePHI[ipar].push_back(
                sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi)));
            m_err_eTHETA[ipar].push_back(
                sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta)));
            m_err_eQOP[ipar].push_back(
                sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)));
            m_err_eT[ipar].push_back(
                sqrt(covariance(Acts::eBoundTime, Acts::eBoundTime)));

            // track parameters pull
            m_pull_eLOC0[ipar].push_back(
                (parameters[Acts::eBoundLoc0] - truthLOC0) /
                sqrt(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)));
            m_pull_eLOC1[ipar].push_back(
                (parameters[Acts::eBoundLoc1] - truthLOC1) /
                sqrt(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)));
            m_pull_ePHI[ipar].push_back(
                resPhi / sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi)));
            m_pull_eTHETA[ipar].push_back(
                (parameters[Acts::eBoundTheta] - truthTHETA) /
                sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta)));
            m_pull_eQOP[ipar].push_back(
                (parameters[Acts::eBoundQOverP] - truthQOP) /
                sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)));
            double sigmaTime =
                sqrt(covariance(Acts::eBoundTime, Acts::eBoundTime));
            m_pull_eT[ipar].push_back(
                sigmaTime == 0.0
                    ? std::numeric_limits<double>::quiet_NaN()
                    : (parameters[Acts::eBoundTime] - truthTIME) / sigmaTime);

            // further track parameter info
            Acts::FreeVector freeParams =
                Acts::detail::transformBoundToFreeParameters(surface, gctx,
                                                             parameters);
            m_x[ipar].push_back(freeParams[Acts::eFreePos0]);
            m_y[ipar].push_back(freeParams[Acts::eFreePos1]);
            m_z[ipar].push_back(freeParams[Acts::eFreePos2]);
            auto p = std::abs(1 / freeParams[Acts::eFreeQOverP]);
            m_px[ipar].push_back(p * freeParams[Acts::eFreeDir0]);
            m_py[ipar].push_back(p * freeParams[Acts::eFreeDir1]);
            m_pz[ipar].push_back(p * freeParams[Acts::eFreeDir2]);
            m_pT[ipar].push_back(p * std::hypot(freeParams[Acts::eFreeDir0],
                                                freeParams[Acts::eFreeDir1]));
            m_eta[ipar].push_back(Acts::VectorHelpers::eta(
                freeParams.segment<3>(Acts::eFreeDir0)));
          } else {
            if (ipar == 0) {
              // push default values if no track parameters
              m_res_x_hit.push_back(-99.);
              m_res_y_hit.push_back(-99.);
              m_err_x_hit.push_back(-99.);
              m_err_y_hit.push_back(-99.);
              m_pull_x_hit.push_back(-99.);
              m_pull_y_hit.push_back(-99.);
              m_dim_hit.push_back(-99.);
            }
            // push default values if no track parameters
            m_eLOC0[ipar].push_back(-99.);
            m_eLOC1[ipar].push_back(-99.);
            m_ePHI[ipar].push_back(-99.);
            m_eTHETA[ipar].push_back(-99.);
            m_eQOP[ipar].push_back(-99.);
            m_eT[ipar].push_back(-99.);
            m_res_eLOC0[ipar].push_back(-99.);
            m_res_eLOC1[ipar].push_back(-99.);
            m_res_ePHI[ipar].push_back(-99.);
            m_res_eTHETA[ipar].push_back(-99.);
            m_res_eQOP[ipar].push_back(-99.);
            m_res_eT[ipar].push_back(-99.);
            m_err_eLOC0[ipar].push_back(-99);
            m_err_eLOC1[ipar].push_back(-99);
            m_err_ePHI[ipar].push_back(-99);
            m_err_eTHETA[ipar].push_back(-99);
            m_err_eQOP[ipar].push_back(-99);
            m_err_eT[ipar].push_back(-99);
            m_pull_eLOC0[ipar].push_back(-99.);
            m_pull_eLOC1[ipar].push_back(-99.);
            m_pull_ePHI[ipar].push_back(-99.);
            m_pull_eTHETA[ipar].push_back(-99.);
            m_pull_eQOP[ipar].push_back(-99.);
            m_pull_eT[ipar].push_back(-99.);
            m_x[ipar].push_back(-99.);
            m_y[ipar].push_back(-99.);
            m_z[ipar].push_back(-99.);
            m_px[ipar].push_back(-99.);
            m_py[ipar].push_back(-99.);
            m_pz[ipar].push_back(-99.);
            m_pT[ipar].push_back(-99.);
            m_eta[ipar].push_back(-99.);
          }
          // fill the track parameters status
          m_hasParams[ipar].push_back(hasParams[ipar]);
        }

        // fill the chi2
        m_chi2.push_back(state.chi2());

        return true;
      });  // all states

      // fill the variables for one track to tree
      m_outputTree->Fill();

      // now reset
      m_t_x.clear();
      m_t_y.clear();
      m_t_z.clear();
      m_t_r.clear();
      m_t_dx.clear();
      m_t_dy.clear();
      m_t_dz.clear();
      m_t_eLOC0.clear();
      m_t_eLOC1.clear();
      m_t_ePHI.clear();
      m_t_eTHETA.clear();
      m_t_eQOP.clear();
      m_t_eT.clear();

      m_volumeID.clear();
      m_layerID.clear();
      m_moduleID.clear();
      m_pathLength.clear();
      m_lx_hit.clear();
      m_ly_hit.clear();
      m_x_hit.clear();
      m_y_hit.clear();
      m_z_hit.clear();
      m_res_x_hit.clear();
      m_res_y_hit.clear();
      m_err_x_hit.clear();
      m_err_y_hit.clear();
      m_pull_x_hit.clear();
      m_pull_y_hit.clear();
      m_dim_hit.clear();

      for (unsigned int ipar = 0; ipar < 3; ++ipar) {
        m_hasParams[ipar].clear();
        m_eLOC0[ipar].clear();
        m_eLOC1[ipar].clear();
        m_ePHI[ipar].clear();
        m_eTHETA[ipar].clear();
        m_eQOP[ipar].clear();
        m_eT[ipar].clear();
        m_res_eLOC0[ipar].clear();
        m_res_eLOC1[ipar].clear();
        m_res_ePHI[ipar].clear();
        m_res_eTHETA[ipar].clear();
        m_res_eQOP[ipar].clear();
        m_res_eT[ipar].clear();
        m_err_eLOC0[ipar].clear();
        m_err_eLOC1[ipar].clear();
        m_err_ePHI[ipar].clear();
        m_err_eTHETA[ipar].clear();
        m_err_eQOP[ipar].clear();
        m_err_eT[ipar].clear();
        m_pull_eLOC0[ipar].clear();
        m_pull_eLOC1[ipar].clear();
        m_pull_ePHI[ipar].clear();
        m_pull_eTHETA[ipar].clear();
        m_pull_eQOP[ipar].clear();
        m_pull_eT[ipar].clear();
        m_x[ipar].clear();
        m_y[ipar].clear();
        m_z[ipar].clear();
        m_px[ipar].clear();
        m_py[ipar].clear();
        m_pz[ipar].clear();
        m_eta[ipar].clear();
        m_pT[ipar].clear();
      }

      m_chi2.clear();
    }  // all subtrajectories
  }    // all trajectories

  return ProcessCode::SUCCESS;
}

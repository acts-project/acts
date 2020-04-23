// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Io/Root/RootTrajectoryWriter.hpp"

#include <TFile.h>
#include <TTree.h>
#include <ios>
#include <stdexcept>

#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Helpers.hpp"

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;
using Measurement = Acts::Measurement<FW::SimSourceLink, Acts::ParDef::eLOC_0,
                                      Acts::ParDef::eLOC_1>;

FW::RootTrajectoryWriter::RootTrajectoryWriter(
    const FW::RootTrajectoryWriter::Config& cfg, Acts::Logging::Level lvl)
    : WriterT(cfg.inputTrajectories, "RootTrajectoryWriter", lvl),
      m_cfg(cfg),
      m_outputFile(cfg.rootFile) {
  // An input collection name and tree name must be specified
  if (m_cfg.inputTrajectories.empty()) {
    throw std::invalid_argument("Missing input trajectory collection");
  } else if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particle collection");
  } else if (cfg.outputFilename.empty()) {
    throw std::invalid_argument("Missing output filename");
  } else if (m_cfg.outputTreename.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // Setup ROOT I/O
  if (m_outputFile == nullptr) {
    auto path = joinPaths(m_cfg.outputDir, m_cfg.outputFilename);
    m_outputFile = TFile::Open(path.c_str(), m_cfg.fileMode.c_str());
    if (m_outputFile == nullptr) {
      throw std::ios_base::failure("Could not open '" + path);
    }
  }
  m_outputFile->cd();
  m_outputTree =
      new TTree(m_cfg.outputTreename.c_str(), m_cfg.outputTreename.c_str());
  if (m_outputTree == nullptr)
    throw std::bad_alloc();
  else {
    // I/O parameters
    m_outputTree->Branch("event_nr", &m_eventNr);
    m_outputTree->Branch("traj_nr", &m_trajNr);
    m_outputTree->Branch("t_barcode", &m_t_barcode, "t_barcode/l");
    m_outputTree->Branch("t_charge", &m_t_charge);
    m_outputTree->Branch("t_time", &m_t_time);
    m_outputTree->Branch("t_vx", &m_t_vx);
    m_outputTree->Branch("t_vy", &m_t_vy);
    m_outputTree->Branch("t_vz", &m_t_vz);
    m_outputTree->Branch("t_px", &m_t_px);
    m_outputTree->Branch("t_py", &m_t_py);
    m_outputTree->Branch("t_pz", &m_t_pz);
    m_outputTree->Branch("t_theta", &m_t_theta);
    m_outputTree->Branch("t_phi", &m_t_phi);
    m_outputTree->Branch("t_eta", &m_t_eta);
    m_outputTree->Branch("t_pT", &m_t_pT);

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

    m_outputTree->Branch("hasFittedParams", &m_hasFittedParams);
    m_outputTree->Branch("eLOC0_fit", &m_eLOC0_fit);
    m_outputTree->Branch("eLOC1_fit", &m_eLOC1_fit);
    m_outputTree->Branch("ePHI_fit", &m_ePHI_fit);
    m_outputTree->Branch("eTHETA_fit", &m_eTHETA_fit);
    m_outputTree->Branch("eQOP_fit", &m_eQOP_fit);
    m_outputTree->Branch("eT_fit", &m_eT_fit);
    m_outputTree->Branch("err_eLOC0_fit", &m_err_eLOC0_fit);
    m_outputTree->Branch("err_eLOC1_fit", &m_err_eLOC1_fit);
    m_outputTree->Branch("err_ePHI_fit", &m_err_ePHI_fit);
    m_outputTree->Branch("err_eTHETA_fit", &m_err_eTHETA_fit);
    m_outputTree->Branch("err_eQOP_fit", &m_err_eQOP_fit);
    m_outputTree->Branch("err_eT_fit", &m_err_eT_fit);

    m_outputTree->Branch("nPredicted", &m_nPredicted);
    m_outputTree->Branch("predicted", &m_prt);
    m_outputTree->Branch("eLOC0_prt", &m_eLOC0_prt);
    m_outputTree->Branch("eLOC1_prt", &m_eLOC1_prt);
    m_outputTree->Branch("ePHI_prt", &m_ePHI_prt);
    m_outputTree->Branch("eTHETA_prt", &m_eTHETA_prt);
    m_outputTree->Branch("eQOP_prt", &m_eQOP_prt);
    m_outputTree->Branch("eT_prt", &m_eT_prt);
    m_outputTree->Branch("res_eLOC0_prt", &m_res_eLOC0_prt);
    m_outputTree->Branch("res_eLOC1_prt", &m_res_eLOC1_prt);
    m_outputTree->Branch("res_ePHI_prt", &m_res_ePHI_prt);
    m_outputTree->Branch("res_eTHETA_prt", &m_res_eTHETA_prt);
    m_outputTree->Branch("res_eQOP_prt", &m_res_eQOP_prt);
    m_outputTree->Branch("res_eT_prt", &m_res_eT_prt);
    m_outputTree->Branch("err_eLOC0_prt", &m_err_eLOC0_prt);
    m_outputTree->Branch("err_eLOC1_prt", &m_err_eLOC1_prt);
    m_outputTree->Branch("err_ePHI_prt", &m_err_ePHI_prt);
    m_outputTree->Branch("err_eTHETA_prt", &m_err_eTHETA_prt);
    m_outputTree->Branch("err_eQOP_prt", &m_err_eQOP_prt);
    m_outputTree->Branch("err_eT_prt", &m_err_eT_prt);
    m_outputTree->Branch("pull_eLOC0_prt", &m_pull_eLOC0_prt);
    m_outputTree->Branch("pull_eLOC1_prt", &m_pull_eLOC1_prt);
    m_outputTree->Branch("pull_ePHI_prt", &m_pull_ePHI_prt);
    m_outputTree->Branch("pull_eTHETA_prt", &m_pull_eTHETA_prt);
    m_outputTree->Branch("pull_eQOP_prt", &m_pull_eQOP_prt);
    m_outputTree->Branch("pull_eT_prt", &m_pull_eT_prt);
    m_outputTree->Branch("g_x_prt", &m_x_prt);
    m_outputTree->Branch("g_y_prt", &m_y_prt);
    m_outputTree->Branch("g_z_prt", &m_z_prt);
    m_outputTree->Branch("px_prt", &m_px_prt);
    m_outputTree->Branch("py_prt", &m_py_prt);
    m_outputTree->Branch("pz_prt", &m_pz_prt);
    m_outputTree->Branch("eta_prt", &m_eta_prt);
    m_outputTree->Branch("pT_prt", &m_pT_prt);

    m_outputTree->Branch("nFiltered", &m_nFiltered);
    m_outputTree->Branch("filtered", &m_flt);
    m_outputTree->Branch("eLOC0_flt", &m_eLOC0_flt);
    m_outputTree->Branch("eLOC1_flt", &m_eLOC1_flt);
    m_outputTree->Branch("ePHI_flt", &m_ePHI_flt);
    m_outputTree->Branch("eTHETA_flt", &m_eTHETA_flt);
    m_outputTree->Branch("eQOP_flt", &m_eQOP_flt);
    m_outputTree->Branch("eT_flt", &m_eT_flt);
    m_outputTree->Branch("res_eLOC0_flt", &m_res_eLOC0_flt);
    m_outputTree->Branch("res_eLOC1_flt", &m_res_eLOC1_flt);
    m_outputTree->Branch("res_ePHI_flt", &m_res_ePHI_flt);
    m_outputTree->Branch("res_eTHETA_flt", &m_res_eTHETA_flt);
    m_outputTree->Branch("res_eQOP_flt", &m_res_eQOP_flt);
    m_outputTree->Branch("res_eT_flt", &m_res_eT_flt);
    m_outputTree->Branch("err_eLOC0_flt", &m_err_eLOC0_flt);
    m_outputTree->Branch("err_eLOC1_flt", &m_err_eLOC1_flt);
    m_outputTree->Branch("err_ePHI_flt", &m_err_ePHI_flt);
    m_outputTree->Branch("err_eTHETA_flt", &m_err_eTHETA_flt);
    m_outputTree->Branch("err_eQOP_flt", &m_err_eQOP_flt);
    m_outputTree->Branch("err_eT_flt", &m_err_eT_flt);
    m_outputTree->Branch("pull_eLOC0_flt", &m_pull_eLOC0_flt);
    m_outputTree->Branch("pull_eLOC1_flt", &m_pull_eLOC1_flt);
    m_outputTree->Branch("pull_ePHI_flt", &m_pull_ePHI_flt);
    m_outputTree->Branch("pull_eTHETA_flt", &m_pull_eTHETA_flt);
    m_outputTree->Branch("pull_eQOP_flt", &m_pull_eQOP_flt);
    m_outputTree->Branch("pull_eT_flt", &m_pull_eT_flt);
    m_outputTree->Branch("g_x_flt", &m_x_flt);
    m_outputTree->Branch("g_y_flt", &m_y_flt);
    m_outputTree->Branch("g_z_flt", &m_z_flt);
    m_outputTree->Branch("px_flt", &m_px_flt);
    m_outputTree->Branch("py_flt", &m_py_flt);
    m_outputTree->Branch("pz_flt", &m_pz_flt);
    m_outputTree->Branch("eta_flt", &m_eta_flt);
    m_outputTree->Branch("pT_flt", &m_pT_flt);
    m_outputTree->Branch("chi2", &m_chi2);

    m_outputTree->Branch("nSmoothed", &m_nSmoothed);
    m_outputTree->Branch("smoothed", &m_smt);
    m_outputTree->Branch("eLOC0_smt", &m_eLOC0_smt);
    m_outputTree->Branch("eLOC1_smt", &m_eLOC1_smt);
    m_outputTree->Branch("ePHI_smt", &m_ePHI_smt);
    m_outputTree->Branch("eTHETA_smt", &m_eTHETA_smt);
    m_outputTree->Branch("eQOP_smt", &m_eQOP_smt);
    m_outputTree->Branch("eT_smt", &m_eT_smt);
    m_outputTree->Branch("res_eLOC0_smt", &m_res_eLOC0_smt);
    m_outputTree->Branch("res_eLOC1_smt", &m_res_eLOC1_smt);
    m_outputTree->Branch("res_ePHI_smt", &m_res_ePHI_smt);
    m_outputTree->Branch("res_eTHETA_smt", &m_res_eTHETA_smt);
    m_outputTree->Branch("res_eQOP_smt", &m_res_eQOP_smt);
    m_outputTree->Branch("res_eT_smt", &m_res_eT_smt);
    m_outputTree->Branch("err_eLOC0_smt", &m_err_eLOC0_smt);
    m_outputTree->Branch("err_eLOC1_smt", &m_err_eLOC1_smt);
    m_outputTree->Branch("err_ePHI_smt", &m_err_ePHI_smt);
    m_outputTree->Branch("err_eTHETA_smt", &m_err_eTHETA_smt);
    m_outputTree->Branch("err_eQOP_smt", &m_err_eQOP_smt);
    m_outputTree->Branch("err_eT_smt", &m_err_eT_smt);
    m_outputTree->Branch("pull_eLOC0_smt", &m_pull_eLOC0_smt);
    m_outputTree->Branch("pull_eLOC1_smt", &m_pull_eLOC1_smt);
    m_outputTree->Branch("pull_ePHI_smt", &m_pull_ePHI_smt);
    m_outputTree->Branch("pull_eTHETA_smt", &m_pull_eTHETA_smt);
    m_outputTree->Branch("pull_eQOP_smt", &m_pull_eQOP_smt);
    m_outputTree->Branch("pull_eT_smt", &m_pull_eT_smt);
    m_outputTree->Branch("g_x_smt", &m_x_smt);
    m_outputTree->Branch("g_y_smt", &m_y_smt);
    m_outputTree->Branch("g_z_smt", &m_z_smt);
    m_outputTree->Branch("px_smt", &m_px_smt);
    m_outputTree->Branch("py_smt", &m_py_smt);
    m_outputTree->Branch("pz_smt", &m_pz_smt);
    m_outputTree->Branch("eta_smt", &m_eta_smt);
    m_outputTree->Branch("pT_smt", &m_pT_smt);
  }
}

FW::RootTrajectoryWriter::~RootTrajectoryWriter() {
  if (m_outputFile) {
    m_outputFile->Close();
  }
}

FW::ProcessCode FW::RootTrajectoryWriter::endRun() {
  if (m_outputFile) {
    m_outputFile->cd();
    m_outputTree->Write();
    ACTS_INFO("Write trajectories to tree '"
              << m_cfg.outputTreename << "' in '"
              << joinPaths(m_cfg.outputDir, m_cfg.outputFilename) << "'");
  }
  return ProcessCode::SUCCESS;
}

FW::ProcessCode FW::RootTrajectoryWriter::writeT(
    const AlgorithmContext& ctx, const TrajectoryContainer& trajectories) {
  if (m_outputFile == nullptr)
    return ProcessCode::SUCCESS;

  auto& gctx = ctx.geoContext;

  // read truth particles from input collection
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Get the event number
  m_eventNr = ctx.eventNumber;

  // Loop over the trajectories
  int iTraj = 0;
  for (const auto& traj : trajectories) {
    /// Collect the information
    m_trajNr = iTraj;

    // Collect number of trackstates with measurements
    m_nMeasurements = traj.numMeasurements();

    // No entry for the track without measurements in the tree
    if (m_nMeasurements == 0) {
      continue;
    }

    // Collect number of all trackstates
    m_nStates = traj.numStates();

    // Get the majority truth particle to this track
    const auto particleHitCount = traj.identifyMajorityParticle();
    if (not particleHitCount.empty()) {
      // Get the barcode of the majority truth particle
      m_t_barcode = particleHitCount.front().particleId.value();
      // Find the truth particle via the barcode
      auto ip = particles.find(m_t_barcode);
      if (ip != particles.end()) {
        const auto& particle = *ip;
        ACTS_DEBUG("Find the truth particle with barcode = " << m_t_barcode);
        // Get the truth particle info at vertex
        const auto p = particle.absMomentum();
        m_t_charge = particle.charge();
        m_t_time = particle.time();
        m_t_vx = particle.position().x();
        m_t_vy = particle.position().y();
        m_t_vz = particle.position().z();
        m_t_px = p * particle.unitDirection().x();
        m_t_py = p * particle.unitDirection().y();
        m_t_pz = p * particle.unitDirection().z();
        m_t_theta = theta(particle.unitDirection());
        m_t_phi = phi(particle.unitDirection());
        m_t_eta = eta(particle.unitDirection());
        m_t_pT = p * perp(particle.unitDirection());
      } else {
        ACTS_WARNING("Truth particle with barcode = " << m_t_barcode
                                                      << " not found!");
      }
    }

    // Get the fitted track parameter
    m_hasFittedParams = false;
    if (traj.hasTrackParameters()) {
      m_hasFittedParams = true;
      const auto& boundParam = traj.trackParameters();
      const auto& parameter = boundParam.parameters();
      const auto& covariance = *boundParam.covariance();
      m_eLOC0_fit = parameter[Acts::ParDef::eLOC_0];
      m_eLOC1_fit = parameter[Acts::ParDef::eLOC_1];
      m_ePHI_fit = parameter[Acts::ParDef::ePHI];
      m_eTHETA_fit = parameter[Acts::ParDef::eTHETA];
      m_eQOP_fit = parameter[Acts::ParDef::eQOP];
      m_eT_fit = parameter[Acts::ParDef::eT];
      m_err_eLOC0_fit =
          sqrt(covariance(Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_0));
      m_err_eLOC1_fit =
          sqrt(covariance(Acts::ParDef::eLOC_1, Acts::ParDef::eLOC_1));
      m_err_ePHI_fit = sqrt(covariance(Acts::ParDef::ePHI, Acts::ParDef::ePHI));
      m_err_eTHETA_fit =
          sqrt(covariance(Acts::ParDef::eTHETA, Acts::ParDef::eTHETA));
      m_err_eQOP_fit = sqrt(covariance(Acts::ParDef::eQOP, Acts::ParDef::eQOP));
      m_err_eT_fit = sqrt(covariance(Acts::ParDef::eT, Acts::ParDef::eT));
    }

    // Get the fitted trajectory
    const auto& [trackTip, mj] = traj.trajectory();

    // Get the trackStates on the trajectory
    m_nPredicted = 0;
    m_nFiltered = 0;
    m_nSmoothed = 0;
    mj.visitBackwards(trackTip, [&](const auto& state) {
      // we only fill the track states with non-outlier measurement
      auto typeFlags = state.typeFlags();
      if (not typeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
        return true;
      }

      // get the geometry ID
      auto geoID = state.referenceSurface().geoID();
      m_volumeID.push_back(geoID.volume());
      m_layerID.push_back(geoID.layer());
      m_moduleID.push_back(geoID.sensitive());

      auto meas = std::get<Measurement>(*state.uncalibrated());

      // get local position
      Acts::Vector2D local(meas.parameters()[Acts::ParDef::eLOC_0],
                           meas.parameters()[Acts::ParDef::eLOC_1]);
      // get global position
      Acts::Vector3D global(0, 0, 0);
      Acts::Vector3D mom(1, 1, 1);
      meas.referenceSurface().localToGlobal(ctx.geoContext, local, mom, global);

      // get measurement covariance
      auto cov = meas.covariance();
      // float resX = sqrt(cov(Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_0));
      // float resY = sqrt(cov(Acts::ParDef::eLOC_1, Acts::ParDef::eLOC_1));

      // push the measurement info
      m_lx_hit.push_back(local.x());
      m_ly_hit.push_back(local.y());
      m_x_hit.push_back(global.x());
      m_y_hit.push_back(global.y());
      m_z_hit.push_back(global.z());

      // get the truth hit corresponding to this trackState
      const auto& truthHit = state.uncalibrated().truthHit();
      // get local truth position
      Acts::Vector2D truthlocal;
      meas.referenceSurface().globalToLocal(
          gctx, truthHit.position(), truthHit.unitDirection(), truthlocal);

      // push the truth hit info
      m_t_x.push_back(truthHit.position().x());
      m_t_y.push_back(truthHit.position().y());
      m_t_z.push_back(truthHit.position().z());
      m_t_r.push_back(perp(truthHit.position()));
      m_t_dx.push_back(truthHit.unitDirection().x());
      m_t_dy.push_back(truthHit.unitDirection().y());
      m_t_dz.push_back(truthHit.unitDirection().z());

      // get the truth track parameter at this track State
      float truthLOC0 = 0, truthLOC1 = 0, truthPHI = 0, truthTHETA = 0,
            truthQOP = 0, truthTIME = 0;
      truthLOC0 = truthlocal.x();
      truthLOC1 = truthlocal.y();
      truthPHI = phi(truthHit.unitDirection());
      truthTHETA = theta(truthHit.unitDirection());
      truthQOP =
          m_t_charge / truthHit.momentum4Before().template head<3>().norm();
      truthTIME = truthHit.time();

      // push the truth track parameter at this track State
      m_t_eLOC0.push_back(truthLOC0);
      m_t_eLOC1.push_back(truthLOC1);
      m_t_ePHI.push_back(truthPHI);
      m_t_eTHETA.push_back(truthTHETA);
      m_t_eQOP.push_back(truthQOP);
      m_t_eT.push_back(truthTIME);

      // get the predicted parameter
      bool predicted = false;
      if (state.hasPredicted()) {
        predicted = true;
        m_nPredicted++;
        Acts::BoundParameters parameter(
            gctx, state.predictedCovariance(), state.predicted(),
            state.referenceSurface().getSharedPtr());
        auto covariance = state.predictedCovariance();
        // local hit residual info
        auto H = meas.projector();
        auto resCov = cov + H * covariance * H.transpose();
        auto residual = meas.residual(parameter);
        m_res_x_hit.push_back(residual(Acts::ParDef::eLOC_0));
        m_res_y_hit.push_back(residual(Acts::ParDef::eLOC_1));
        m_err_x_hit.push_back(
            sqrt(resCov(Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_0)));
        m_err_y_hit.push_back(
            sqrt(resCov(Acts::ParDef::eLOC_1, Acts::ParDef::eLOC_1)));
        m_pull_x_hit.push_back(
            residual(Acts::ParDef::eLOC_0) /
            sqrt(resCov(Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_0)));
        m_pull_y_hit.push_back(
            residual(Acts::ParDef::eLOC_1) /
            sqrt(resCov(Acts::ParDef::eLOC_1, Acts::ParDef::eLOC_1)));
        m_dim_hit.push_back(state.calibratedSize());

        // predicted parameter
        m_eLOC0_prt.push_back(parameter.parameters()[Acts::ParDef::eLOC_0]);
        m_eLOC1_prt.push_back(parameter.parameters()[Acts::ParDef::eLOC_1]);
        m_ePHI_prt.push_back(parameter.parameters()[Acts::ParDef::ePHI]);
        m_eTHETA_prt.push_back(parameter.parameters()[Acts::ParDef::eTHETA]);
        m_eQOP_prt.push_back(parameter.parameters()[Acts::ParDef::eQOP]);
        m_eT_prt.push_back(parameter.parameters()[Acts::ParDef::eT]);

        // predicted residual
        m_res_eLOC0_prt.push_back(parameter.parameters()[Acts::ParDef::eLOC_0] -
                                  truthLOC0);
        m_res_eLOC1_prt.push_back(parameter.parameters()[Acts::ParDef::eLOC_1] -
                                  truthLOC1);
        m_res_ePHI_prt.push_back(parameter.parameters()[Acts::ParDef::ePHI] -
                                 truthPHI);
        m_res_eTHETA_prt.push_back(
            parameter.parameters()[Acts::ParDef::eTHETA] - truthTHETA);
        m_res_eQOP_prt.push_back(parameter.parameters()[Acts::ParDef::eQOP] -
                                 truthQOP);
        m_res_eT_prt.push_back(parameter.parameters()[Acts::ParDef::eT] -
                               truthTIME);

        // predicted parameter error
        m_err_eLOC0_prt.push_back(
            sqrt(covariance(Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_0)));
        m_err_eLOC1_prt.push_back(
            sqrt(covariance(Acts::ParDef::eLOC_1, Acts::ParDef::eLOC_1)));
        m_err_ePHI_prt.push_back(
            sqrt(covariance(Acts::ParDef::ePHI, Acts::ParDef::ePHI)));
        m_err_eTHETA_prt.push_back(
            sqrt(covariance(Acts::ParDef::eTHETA, Acts::ParDef::eTHETA)));
        m_err_eQOP_prt.push_back(
            sqrt(covariance(Acts::ParDef::eQOP, Acts::ParDef::eQOP)));
        m_err_eT_prt.push_back(
            sqrt(covariance(Acts::ParDef::eT, Acts::ParDef::eT)));

        // predicted parameter pull
        m_pull_eLOC0_prt.push_back(
            (parameter.parameters()[Acts::ParDef::eLOC_0] - truthLOC0) /
            sqrt(covariance(Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_0)));
        m_pull_eLOC1_prt.push_back(
            (parameter.parameters()[Acts::ParDef::eLOC_1] - truthLOC1) /
            sqrt(covariance(Acts::ParDef::eLOC_1, Acts::ParDef::eLOC_1)));
        m_pull_ePHI_prt.push_back(
            (parameter.parameters()[Acts::ParDef::ePHI] - truthPHI) /
            sqrt(covariance(Acts::ParDef::ePHI, Acts::ParDef::ePHI)));
        m_pull_eTHETA_prt.push_back(
            (parameter.parameters()[Acts::ParDef::eTHETA] - truthTHETA) /
            sqrt(covariance(Acts::ParDef::eTHETA, Acts::ParDef::eTHETA)));
        m_pull_eQOP_prt.push_back(
            (parameter.parameters()[Acts::ParDef::eQOP] - truthQOP) /
            sqrt(covariance(Acts::ParDef::eQOP, Acts::ParDef::eQOP)));
        m_pull_eT_prt.push_back(
            (parameter.parameters()[Acts::ParDef::eT] - truthTIME) /
            sqrt(covariance(Acts::ParDef::eT, Acts::ParDef::eT)));

        // further predicted parameter info
        m_x_prt.push_back(parameter.position().x());
        m_y_prt.push_back(parameter.position().y());
        m_z_prt.push_back(parameter.position().z());
        m_px_prt.push_back(parameter.momentum().x());
        m_py_prt.push_back(parameter.momentum().y());
        m_pz_prt.push_back(parameter.momentum().z());
        m_pT_prt.push_back(parameter.pT());
        m_eta_prt.push_back(eta(parameter.position()));
      } else {
        // push default values if no predicted parameter
        m_res_x_hit.push_back(-99.);
        m_res_y_hit.push_back(-99.);
        m_err_x_hit.push_back(-99.);
        m_err_y_hit.push_back(-99.);
        m_pull_x_hit.push_back(-99.);
        m_pull_y_hit.push_back(-99.);
        m_dim_hit.push_back(-99.);
        m_eLOC0_prt.push_back(-99.);
        m_eLOC1_prt.push_back(-99.);
        m_ePHI_prt.push_back(-99.);
        m_eTHETA_prt.push_back(-99.);
        m_eQOP_prt.push_back(-99.);
        m_eT_prt.push_back(-99.);
        m_res_eLOC0_prt.push_back(-99.);
        m_res_eLOC1_prt.push_back(-99.);
        m_res_ePHI_prt.push_back(-99.);
        m_res_eTHETA_prt.push_back(-99.);
        m_res_eQOP_prt.push_back(-99.);
        m_res_eT_prt.push_back(-99.);
        m_err_eLOC0_prt.push_back(-99);
        m_err_eLOC1_prt.push_back(-99);
        m_err_ePHI_prt.push_back(-99);
        m_err_eTHETA_prt.push_back(-99);
        m_err_eQOP_prt.push_back(-99);
        m_err_eT_prt.push_back(-99);
        m_pull_eLOC0_prt.push_back(-99.);
        m_pull_eLOC1_prt.push_back(-99.);
        m_pull_ePHI_prt.push_back(-99.);
        m_pull_eTHETA_prt.push_back(-99.);
        m_pull_eQOP_prt.push_back(-99.);
        m_pull_eT_prt.push_back(-99.);
        m_x_prt.push_back(-99.);
        m_y_prt.push_back(-99.);
        m_z_prt.push_back(-99.);
        m_px_prt.push_back(-99.);
        m_py_prt.push_back(-99.);
        m_pz_prt.push_back(-99.);
        m_pT_prt.push_back(-99.);
        m_eta_prt.push_back(-99.);
      }

      // get the filtered parameter
      bool filtered = false;
      if (state.hasFiltered()) {
        filtered = true;
        m_nFiltered++;
        Acts::BoundParameters parameter(
            gctx, state.filteredCovariance(), state.filtered(),
            state.referenceSurface().getSharedPtr());
        auto covariance = state.filteredCovariance();
        // filtered parameter
        m_eLOC0_flt.push_back(parameter.parameters()[Acts::ParDef::eLOC_0]);
        m_eLOC1_flt.push_back(parameter.parameters()[Acts::ParDef::eLOC_1]);
        m_ePHI_flt.push_back(parameter.parameters()[Acts::ParDef::ePHI]);
        m_eTHETA_flt.push_back(parameter.parameters()[Acts::ParDef::eTHETA]);
        m_eQOP_flt.push_back(parameter.parameters()[Acts::ParDef::eQOP]);
        m_eT_flt.push_back(parameter.parameters()[Acts::ParDef::eT]);

        // filtered residual
        m_res_eLOC0_flt.push_back(parameter.parameters()[Acts::ParDef::eLOC_0] -
                                  truthLOC0);
        m_res_eLOC1_flt.push_back(parameter.parameters()[Acts::ParDef::eLOC_1] -
                                  truthLOC1);
        m_res_ePHI_flt.push_back(parameter.parameters()[Acts::ParDef::ePHI] -
                                 truthPHI);
        m_res_eTHETA_flt.push_back(
            parameter.parameters()[Acts::ParDef::eTHETA] - truthTHETA);
        m_res_eQOP_flt.push_back(parameter.parameters()[Acts::ParDef::eQOP] -
                                 truthQOP);
        m_res_eT_flt.push_back(parameter.parameters()[Acts::ParDef::eT] -
                               truthTIME);

        // filtered parameter error
        m_err_eLOC0_flt.push_back(
            sqrt(covariance(Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_0)));
        m_err_eLOC1_flt.push_back(
            sqrt(covariance(Acts::ParDef::eLOC_1, Acts::ParDef::eLOC_1)));
        m_err_ePHI_flt.push_back(
            sqrt(covariance(Acts::ParDef::ePHI, Acts::ParDef::ePHI)));
        m_err_eTHETA_flt.push_back(
            sqrt(covariance(Acts::ParDef::eTHETA, Acts::ParDef::eTHETA)));
        m_err_eQOP_flt.push_back(
            sqrt(covariance(Acts::ParDef::eQOP, Acts::ParDef::eQOP)));
        m_err_eT_flt.push_back(
            sqrt(covariance(Acts::ParDef::eT, Acts::ParDef::eT)));

        // filtered parameter pull
        m_pull_eLOC0_flt.push_back(
            (parameter.parameters()[Acts::ParDef::eLOC_0] - truthLOC0) /
            sqrt(covariance(Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_0)));
        m_pull_eLOC1_flt.push_back(
            (parameter.parameters()[Acts::ParDef::eLOC_1] - truthLOC1) /
            sqrt(covariance(Acts::ParDef::eLOC_1, Acts::ParDef::eLOC_1)));
        m_pull_ePHI_flt.push_back(
            (parameter.parameters()[Acts::ParDef::ePHI] - truthPHI) /
            sqrt(covariance(Acts::ParDef::ePHI, Acts::ParDef::ePHI)));
        m_pull_eTHETA_flt.push_back(
            (parameter.parameters()[Acts::ParDef::eTHETA] - truthTHETA) /
            sqrt(covariance(Acts::ParDef::eTHETA, Acts::ParDef::eTHETA)));
        m_pull_eQOP_flt.push_back(
            (parameter.parameters()[Acts::ParDef::eQOP] - truthQOP) /
            sqrt(covariance(Acts::ParDef::eQOP, Acts::ParDef::eQOP)));
        m_pull_eT_flt.push_back(
            (parameter.parameters()[Acts::ParDef::eT] - truthTIME) /
            sqrt(covariance(Acts::ParDef::eT, Acts::ParDef::eT)));

        // more filtered parameter info
        m_x_flt.push_back(parameter.position().x());
        m_y_flt.push_back(parameter.position().y());
        m_z_flt.push_back(parameter.position().z());
        m_px_flt.push_back(parameter.momentum().x());
        m_py_flt.push_back(parameter.momentum().y());
        m_pz_flt.push_back(parameter.momentum().z());
        m_pT_flt.push_back(parameter.pT());
        m_eta_flt.push_back(eta(parameter.position()));
        m_chi2.push_back(state.chi2());
      } else {
        // push default values if no filtered parameter
        m_eLOC0_flt.push_back(-99.);
        m_eLOC1_flt.push_back(-99.);
        m_ePHI_flt.push_back(-99.);
        m_eTHETA_flt.push_back(-99.);
        m_eQOP_flt.push_back(-99.);
        m_eT_flt.push_back(-99.);
        m_res_eLOC0_flt.push_back(-99.);
        m_res_eLOC1_flt.push_back(-99.);
        m_res_ePHI_flt.push_back(-99.);
        m_res_eTHETA_flt.push_back(-99.);
        m_res_eQOP_flt.push_back(-99.);
        m_res_eT_flt.push_back(-99.);
        m_err_eLOC0_flt.push_back(-99);
        m_err_eLOC1_flt.push_back(-99);
        m_err_ePHI_flt.push_back(-99);
        m_err_eTHETA_flt.push_back(-99);
        m_err_eQOP_flt.push_back(-99);
        m_err_eT_flt.push_back(-99);
        m_pull_eLOC0_flt.push_back(-99.);
        m_pull_eLOC1_flt.push_back(-99.);
        m_pull_ePHI_flt.push_back(-99.);
        m_pull_eTHETA_flt.push_back(-99.);
        m_pull_eQOP_flt.push_back(-99.);
        m_pull_eT_flt.push_back(-99.);
        m_x_flt.push_back(-99.);
        m_y_flt.push_back(-99.);
        m_z_flt.push_back(-99.);
        m_py_flt.push_back(-99.);
        m_pz_flt.push_back(-99.);
        m_pT_flt.push_back(-99.);
        m_eta_flt.push_back(-99.);
        m_chi2.push_back(-99.0);
      }

      // get the smoothed parameter
      bool smoothed = false;
      if (state.hasSmoothed()) {
        smoothed = true;
        m_nSmoothed++;
        Acts::BoundParameters parameter(
            gctx, state.smoothedCovariance(), state.smoothed(),
            state.referenceSurface().getSharedPtr());
        auto covariance = state.smoothedCovariance();

        // smoothed parameter
        m_eLOC0_smt.push_back(parameter.parameters()[Acts::ParDef::eLOC_0]);
        m_eLOC1_smt.push_back(parameter.parameters()[Acts::ParDef::eLOC_1]);
        m_ePHI_smt.push_back(parameter.parameters()[Acts::ParDef::ePHI]);
        m_eTHETA_smt.push_back(parameter.parameters()[Acts::ParDef::eTHETA]);
        m_eQOP_smt.push_back(parameter.parameters()[Acts::ParDef::eQOP]);
        m_eT_smt.push_back(parameter.parameters()[Acts::ParDef::eT]);

        // smoothed residual
        m_res_eLOC0_smt.push_back(parameter.parameters()[Acts::ParDef::eLOC_0] -
                                  truthLOC0);
        m_res_eLOC1_smt.push_back(parameter.parameters()[Acts::ParDef::eLOC_1] -
                                  truthLOC1);
        m_res_ePHI_smt.push_back(parameter.parameters()[Acts::ParDef::ePHI] -
                                 truthPHI);
        m_res_eTHETA_smt.push_back(
            parameter.parameters()[Acts::ParDef::eTHETA] - truthTHETA);
        m_res_eQOP_smt.push_back(parameter.parameters()[Acts::ParDef::eQOP] -
                                 truthQOP);
        m_res_eT_smt.push_back(parameter.parameters()[Acts::ParDef::eT] -
                               truthTIME);

        // smoothed parameter error
        m_err_eLOC0_smt.push_back(
            sqrt(covariance(Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_0)));
        m_err_eLOC1_smt.push_back(
            sqrt(covariance(Acts::ParDef::eLOC_1, Acts::ParDef::eLOC_1)));
        m_err_ePHI_smt.push_back(
            sqrt(covariance(Acts::ParDef::ePHI, Acts::ParDef::ePHI)));
        m_err_eTHETA_smt.push_back(
            sqrt(covariance(Acts::ParDef::eTHETA, Acts::ParDef::eTHETA)));
        m_err_eQOP_smt.push_back(
            sqrt(covariance(Acts::ParDef::eQOP, Acts::ParDef::eQOP)));
        m_err_eT_smt.push_back(
            sqrt(covariance(Acts::ParDef::eT, Acts::ParDef::eT)));

        // smoothed parameter pull
        m_pull_eLOC0_smt.push_back(
            (parameter.parameters()[Acts::ParDef::eLOC_0] - truthLOC0) /
            sqrt(covariance(Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_0)));
        m_pull_eLOC1_smt.push_back(
            (parameter.parameters()[Acts::ParDef::eLOC_1] - truthLOC1) /
            sqrt(covariance(Acts::ParDef::eLOC_1, Acts::ParDef::eLOC_1)));
        m_pull_ePHI_smt.push_back(
            (parameter.parameters()[Acts::ParDef::ePHI] - truthPHI) /
            sqrt(covariance(Acts::ParDef::ePHI, Acts::ParDef::ePHI)));
        m_pull_eTHETA_smt.push_back(
            (parameter.parameters()[Acts::ParDef::eTHETA] - truthTHETA) /
            sqrt(covariance(Acts::ParDef::eTHETA, Acts::ParDef::eTHETA)));
        m_pull_eQOP_smt.push_back(
            (parameter.parameters()[Acts::ParDef::eQOP] - truthQOP) /
            sqrt(covariance(Acts::ParDef::eQOP, Acts::ParDef::eQOP)));
        m_pull_eT_smt.push_back(
            (parameter.parameters()[Acts::ParDef::eT] - truthTIME) /
            sqrt(covariance(Acts::ParDef::eT, Acts::ParDef::eT)));

        // further smoothed parameter info
        m_x_smt.push_back(parameter.position().x());
        m_y_smt.push_back(parameter.position().y());
        m_z_smt.push_back(parameter.position().z());
        m_px_smt.push_back(parameter.momentum().x());
        m_py_smt.push_back(parameter.momentum().y());
        m_pz_smt.push_back(parameter.momentum().z());
        m_pT_smt.push_back(parameter.pT());
        m_eta_smt.push_back(eta(parameter.position()));
      } else {
        // push default values if no smoothed parameter
        m_eLOC0_smt.push_back(-99.);
        m_eLOC1_smt.push_back(-99.);
        m_ePHI_smt.push_back(-99.);
        m_eTHETA_smt.push_back(-99.);
        m_eQOP_smt.push_back(-99.);
        m_eT_smt.push_back(-99.);
        m_res_eLOC0_smt.push_back(-99.);
        m_res_eLOC1_smt.push_back(-99.);
        m_res_ePHI_smt.push_back(-99.);
        m_res_eTHETA_smt.push_back(-99.);
        m_res_eQOP_smt.push_back(-99.);
        m_res_eT_smt.push_back(-99.);
        m_err_eLOC0_smt.push_back(-99);
        m_err_eLOC1_smt.push_back(-99);
        m_err_ePHI_smt.push_back(-99);
        m_err_eTHETA_smt.push_back(-99);
        m_err_eQOP_smt.push_back(-99);
        m_err_eT_smt.push_back(-99);
        m_pull_eLOC0_smt.push_back(-99.);
        m_pull_eLOC1_smt.push_back(-99.);
        m_pull_ePHI_smt.push_back(-99.);
        m_pull_eTHETA_smt.push_back(-99.);
        m_pull_eQOP_smt.push_back(-99.);
        m_pull_eT_smt.push_back(-99.);
        m_x_smt.push_back(-99.);
        m_y_smt.push_back(-99.);
        m_z_smt.push_back(-99.);
        m_px_smt.push_back(-99.);
        m_py_smt.push_back(-99.);
        m_pz_smt.push_back(-99.);
        m_pT_smt.push_back(-99.);
        m_eta_smt.push_back(-99.);
      }

      m_prt.push_back(predicted);
      m_flt.push_back(filtered);
      m_smt.push_back(smoothed);
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

    m_prt.clear();
    m_eLOC0_prt.clear();
    m_eLOC1_prt.clear();
    m_ePHI_prt.clear();
    m_eTHETA_prt.clear();
    m_eQOP_prt.clear();
    m_eT_prt.clear();
    m_res_eLOC0_prt.clear();
    m_res_eLOC1_prt.clear();
    m_res_ePHI_prt.clear();
    m_res_eTHETA_prt.clear();
    m_res_eQOP_prt.clear();
    m_res_eT_prt.clear();
    m_err_eLOC0_prt.clear();
    m_err_eLOC1_prt.clear();
    m_err_ePHI_prt.clear();
    m_err_eTHETA_prt.clear();
    m_err_eQOP_prt.clear();
    m_err_eT_prt.clear();
    m_pull_eLOC0_prt.clear();
    m_pull_eLOC1_prt.clear();
    m_pull_ePHI_prt.clear();
    m_pull_eTHETA_prt.clear();
    m_pull_eQOP_prt.clear();
    m_pull_eT_prt.clear();
    m_x_prt.clear();
    m_y_prt.clear();
    m_z_prt.clear();
    m_px_prt.clear();
    m_py_prt.clear();
    m_pz_prt.clear();
    m_eta_prt.clear();
    m_pT_prt.clear();

    m_flt.clear();
    m_eLOC0_flt.clear();
    m_eLOC1_flt.clear();
    m_ePHI_flt.clear();
    m_eTHETA_flt.clear();
    m_eQOP_flt.clear();
    m_eT_flt.clear();
    m_res_eLOC0_flt.clear();
    m_res_eLOC1_flt.clear();
    m_res_ePHI_flt.clear();
    m_res_eTHETA_flt.clear();
    m_res_eQOP_flt.clear();
    m_res_eT_flt.clear();
    m_err_eLOC0_flt.clear();
    m_err_eLOC1_flt.clear();
    m_err_ePHI_flt.clear();
    m_err_eTHETA_flt.clear();
    m_err_eQOP_flt.clear();
    m_err_eT_flt.clear();
    m_pull_eLOC0_flt.clear();
    m_pull_eLOC1_flt.clear();
    m_pull_ePHI_flt.clear();
    m_pull_eTHETA_flt.clear();
    m_pull_eQOP_flt.clear();
    m_pull_eT_flt.clear();
    m_x_flt.clear();
    m_y_flt.clear();
    m_z_flt.clear();
    m_px_flt.clear();
    m_py_flt.clear();
    m_pz_flt.clear();
    m_eta_flt.clear();
    m_pT_flt.clear();
    m_chi2.clear();

    m_smt.clear();
    m_eLOC0_smt.clear();
    m_eLOC1_smt.clear();
    m_ePHI_smt.clear();
    m_eTHETA_smt.clear();
    m_eQOP_smt.clear();
    m_eT_smt.clear();
    m_res_eLOC0_smt.clear();
    m_res_eLOC1_smt.clear();
    m_res_ePHI_smt.clear();
    m_res_eTHETA_smt.clear();
    m_res_eQOP_smt.clear();
    m_res_eT_smt.clear();
    m_err_eLOC0_smt.clear();
    m_err_eLOC1_smt.clear();
    m_err_ePHI_smt.clear();
    m_err_eTHETA_smt.clear();
    m_err_eQOP_smt.clear();
    m_err_eT_smt.clear();
    m_pull_eLOC0_smt.clear();
    m_pull_eLOC1_smt.clear();
    m_pull_ePHI_smt.clear();
    m_pull_eTHETA_smt.clear();
    m_pull_eQOP_smt.clear();
    m_pull_eT_smt.clear();
    m_x_smt.clear();
    m_y_smt.clear();
    m_z_smt.clear();
    m_px_smt.clear();
    m_py_smt.clear();
    m_pz_smt.clear();
    m_eta_smt.clear();
    m_pT_smt.clear();

    iTraj++;
  }  // all trajectories

  return ProcessCode::SUCCESS;
}

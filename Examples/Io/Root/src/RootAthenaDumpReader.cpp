// This file is part of the Acts project.
//
// Copyright (C) 2022-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootAthenaDumpReader.hpp"

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Zip.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"

#include <TChain.h>
#include <boost/container/static_vector.hpp>

enum SpacePointType { ePixel = 1, eStrip = 2 };

ActsExamples::RootAthenaDumpReader::RootAthenaDumpReader(
    const ActsExamples::RootAthenaDumpReader::Config& config,
    Acts::Logging::Level level)
    : ActsExamples::IReader(),
      m_cfg(config),
      m_logger(Acts::getDefaultLogger(name(), level)) {
  if (m_cfg.inputfile.empty()) {
    throw std::invalid_argument("Missing input filename");
  }
  if (m_cfg.treename.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  m_inputchain = std::make_shared<TChain>(m_cfg.treename.c_str());

  m_outputPixelSpacePoints.initialize(m_cfg.outputPixelSpacePoints);
  m_outputStripSpacePoints.initialize(m_cfg.outputStripSpacePoints);
  m_outputSpacePoints.initialize(m_cfg.outputSpacePoints);
  m_outputClusters.initialize(m_cfg.outputClusters);

  // Set the branches

  // Set object pointer
  CLhardware = 0;
  CLparticleLink_eventIndex = 0;
  CLparticleLink_barcode = 0;
  CLbarcodesLinked = 0;
  CLparticle_charge = 0;
  CLphis = 0;
  CLetas = 0;
  CLtots = 0;
  CLlocal_cov = 0;
  Part_vParentID = 0;
  Part_vParentBarcode = 0;
  SPtopStripDirection = nullptr;
  SPbottomStripDirection = nullptr;
  SPstripCenterDistance = nullptr;
  SPtopStripCenterPosition = nullptr;
  TRKproperties = 0;
  TRKpattern = 0;
  TRKmeasurementsOnTrack_pixcl_sctcl_index = 0;
  TRKoutliersOnTrack_pixcl_sctcl_index = 0;
  TRKperigee_position = 0;
  TRKperigee_momentum = 0;
  DTTtrajectory_eventindex = 0;
  DTTtrajectory_barcode = 0;
  DTTstTruth_subDetType = 0;
  DTTstTrack_subDetType = 0;
  DTTstCommon_subDetType = 0;

  m_inputchain->SetBranchAddress("run_number", &run_number);
  m_inputchain->SetBranchAddress("event_number", &event_number);
  m_inputchain->SetBranchAddress("nSE", &nSE);
  m_inputchain->SetBranchAddress("SEID", SEID);
  m_inputchain->SetBranchAddress("nCL", &nCL);
  m_inputchain->SetBranchAddress("CLindex", CLindex);
  m_inputchain->SetBranchAddress("CLhardware", &CLhardware);
  m_inputchain->SetBranchAddress("CLx", CLx);
  m_inputchain->SetBranchAddress("CLy", CLy);
  m_inputchain->SetBranchAddress("CLz", CLz);
  m_inputchain->SetBranchAddress("CLbarrel_endcap", CLbarrel_endcap);
  m_inputchain->SetBranchAddress("CLlayer_disk", CLlayer_disk);
  m_inputchain->SetBranchAddress("CLeta_module", CLeta_module);
  m_inputchain->SetBranchAddress("CLphi_module", CLphi_module);
  m_inputchain->SetBranchAddress("CLside", CLside);
  m_inputchain->SetBranchAddress("CLmoduleID", CLmoduleID);
  m_inputchain->SetBranchAddress("CLparticleLink_eventIndex",
                                 &CLparticleLink_eventIndex);
  m_inputchain->SetBranchAddress("CLparticleLink_barcode",
                                 &CLparticleLink_barcode);
  m_inputchain->SetBranchAddress("CLbarcodesLinked", &CLbarcodesLinked);
  m_inputchain->SetBranchAddress("CLparticle_charge", &CLparticle_charge);
  m_inputchain->SetBranchAddress("CLphis", &CLphis);
  m_inputchain->SetBranchAddress("CLetas", &CLetas);
  m_inputchain->SetBranchAddress("CLtots", &CLtots);
  m_inputchain->SetBranchAddress("CLloc_direction1", CLloc_direction1);
  m_inputchain->SetBranchAddress("CLloc_direction2", CLloc_direction2);
  m_inputchain->SetBranchAddress("CLloc_direction3", CLloc_direction3);
  m_inputchain->SetBranchAddress("CLJan_loc_direction1", CLJan_loc_direction1);
  m_inputchain->SetBranchAddress("CLJan_loc_direction2", CLJan_loc_direction2);
  m_inputchain->SetBranchAddress("CLJan_loc_direction3", CLJan_loc_direction3);
  m_inputchain->SetBranchAddress("CLpixel_count", CLpixel_count);
  m_inputchain->SetBranchAddress("CLcharge_count", CLcharge_count);
  m_inputchain->SetBranchAddress("CLloc_eta", CLloc_eta);
  m_inputchain->SetBranchAddress("CLloc_phi", CLloc_phi);
  m_inputchain->SetBranchAddress("CLglob_eta", CLglob_eta);
  m_inputchain->SetBranchAddress("CLglob_phi", CLglob_phi);
  m_inputchain->SetBranchAddress("CLeta_angle", CLeta_angle);
  m_inputchain->SetBranchAddress("CLphi_angle", CLphi_angle);
  m_inputchain->SetBranchAddress("CLnorm_x", CLnorm_x);
  m_inputchain->SetBranchAddress("CLnorm_y", CLnorm_y);
  m_inputchain->SetBranchAddress("CLnorm_z", CLnorm_z);
  m_inputchain->SetBranchAddress("CLlocal_cov", &CLlocal_cov);
  m_inputchain->SetBranchAddress("nPartEVT", &nPartEVT);
  m_inputchain->SetBranchAddress("Part_event_number", Part_event_number);
  m_inputchain->SetBranchAddress("Part_barcode", Part_barcode);
  m_inputchain->SetBranchAddress("Part_px", Part_px);
  m_inputchain->SetBranchAddress("Part_py", Part_py);
  m_inputchain->SetBranchAddress("Part_pz", Part_pz);
  m_inputchain->SetBranchAddress("Part_pt", Part_pt);
  m_inputchain->SetBranchAddress("Part_eta", Part_eta);
  m_inputchain->SetBranchAddress("Part_vx", Part_vx);
  m_inputchain->SetBranchAddress("Part_vy", Part_vy);
  m_inputchain->SetBranchAddress("Part_vz", Part_vz);
  m_inputchain->SetBranchAddress("Part_radius", Part_radius);
  m_inputchain->SetBranchAddress("Part_status", Part_status);
  m_inputchain->SetBranchAddress("Part_charge", Part_charge);
  m_inputchain->SetBranchAddress("Part_pdg_id", Part_pdg_id);
  m_inputchain->SetBranchAddress("Part_passed", Part_passed);
  m_inputchain->SetBranchAddress("Part_vProdNin", Part_vProdNin);
  m_inputchain->SetBranchAddress("Part_vProdNout", Part_vProdNout);
  m_inputchain->SetBranchAddress("Part_vProdStatus", Part_vProdStatus);
  m_inputchain->SetBranchAddress("Part_vProdBarcode", Part_vProdBarcode);
  m_inputchain->SetBranchAddress("Part_vParentID", &Part_vParentID);
  m_inputchain->SetBranchAddress("Part_vParentBarcode", &Part_vParentBarcode);
  m_inputchain->SetBranchAddress("nSP", &nSP);
  m_inputchain->SetBranchAddress("SPindex", SPindex);
  m_inputchain->SetBranchAddress("SPx", SPx);
  m_inputchain->SetBranchAddress("SPy", SPy);
  m_inputchain->SetBranchAddress("SPz", SPz);
  m_inputchain->SetBranchAddress("SPCL1_index", SPCL1_index);
  m_inputchain->SetBranchAddress("SPCL2_index", SPCL2_index);
  m_inputchain->SetBranchAddress("SPisOverlap", SPisOverlap);
  m_inputchain->SetBranchAddress("SPradius", SPradius);
  m_inputchain->SetBranchAddress("SPcovr", SPcovr);
  m_inputchain->SetBranchAddress("SPcovz", SPcovz);
  m_inputchain->SetBranchAddress("SPhl_topstrip", SPhl_topstrip);
  m_inputchain->SetBranchAddress("SPhl_botstrip", SPhl_botstrip);
  m_inputchain->SetBranchAddress("SPtopStripDirection", SPtopStripDirection);
  m_inputchain->SetBranchAddress("SPbottomStripDirection",
                                 SPbottomStripDirection);
  m_inputchain->SetBranchAddress("SPstripCenterDistance",
                                 SPstripCenterDistance);
  m_inputchain->SetBranchAddress("SPtopStripCenterPosition",
                                 SPtopStripCenterPosition);

  m_inputchain->SetBranchAddress("nTRK", &nTRK);
  m_inputchain->SetBranchAddress("TRKindex", TRKindex);
  m_inputchain->SetBranchAddress("TRKtrack_fitter", TRKtrack_fitter);
  m_inputchain->SetBranchAddress("TRKparticle_hypothesis",
                                 TRKparticle_hypothesis);
  m_inputchain->SetBranchAddress("TRKproperties", &TRKproperties);
  m_inputchain->SetBranchAddress("TRKpattern", &TRKpattern);
  m_inputchain->SetBranchAddress("TRKndof", TRKndof);
  m_inputchain->SetBranchAddress("TRKmot", TRKmot);
  m_inputchain->SetBranchAddress("TRKoot", TRKoot);
  m_inputchain->SetBranchAddress("TRKchiSq", TRKchiSq);
  m_inputchain->SetBranchAddress("TRKmeasurementsOnTrack_pixcl_sctcl_index",
                                 &TRKmeasurementsOnTrack_pixcl_sctcl_index);
  m_inputchain->SetBranchAddress("TRKoutliersOnTrack_pixcl_sctcl_index",
                                 &TRKoutliersOnTrack_pixcl_sctcl_index);
  m_inputchain->SetBranchAddress("TRKcharge", TRKcharge);
  m_inputchain->SetBranchAddress("TRKperigee_position", &TRKperigee_position);
  m_inputchain->SetBranchAddress("TRKperigee_momentum", &TRKperigee_momentum);
  m_inputchain->SetBranchAddress("TTCindex", TTCindex);
  m_inputchain->SetBranchAddress("TTCevent_index", TTCevent_index);
  m_inputchain->SetBranchAddress("TTCparticle_link", TTCparticle_link);
  m_inputchain->SetBranchAddress("TTCprobability", TTCprobability);
  m_inputchain->SetBranchAddress("nDTT", &nDTT);
  m_inputchain->SetBranchAddress("DTTindex", DTTindex);
  m_inputchain->SetBranchAddress("DTTsize", DTTsize);
  m_inputchain->SetBranchAddress("DTTtrajectory_eventindex",
                                 &DTTtrajectory_eventindex);
  m_inputchain->SetBranchAddress("DTTtrajectory_barcode",
                                 &DTTtrajectory_barcode);
  m_inputchain->SetBranchAddress("DTTstTruth_subDetType",
                                 &DTTstTruth_subDetType);
  m_inputchain->SetBranchAddress("DTTstTrack_subDetType",
                                 &DTTstTrack_subDetType);
  m_inputchain->SetBranchAddress("DTTstCommon_subDetType",
                                 &DTTstCommon_subDetType);

  m_inputchain->Add(m_cfg.inputfile.c_str());
  ACTS_DEBUG("Adding file " << m_cfg.inputfile << " to tree" << m_cfg.treename);

  m_events = m_inputchain->GetEntries();

  ACTS_DEBUG("End of constructor. In total available events=" << m_events);

}  // constructor

ActsExamples::ProcessCode ActsExamples::RootAthenaDumpReader::read(
    const ActsExamples::AlgorithmContext& ctx) {
  ACTS_DEBUG("Reading event " << ctx.eventNumber);
  auto entry = ctx.eventNumber;
  if (entry >= m_events) {
    ACTS_ERROR("event out of bounds");
    return ProcessCode::ABORT;
  }

  std::lock_guard<std::mutex> lock(m_read_mutex);

  m_inputchain->GetEntry(entry);

  // Loop on clusters (measurements)
  ACTS_DEBUG("Found " << nSP << " space points");
  ACTS_DEBUG("Found " << nCL << " clusters / measurements");

  ClusterContainer clusters;
  clusters.resize(nCL);

  for (int im = 0; im < nCL; im++) {
    int bec = CLbarrel_endcap[im];
    int lydisk = CLlayer_disk[im];
    int etamod = CLeta_module[im];
    int phimod = CLphi_module[im];
    int side = CLside[im];
    // ULong64_t moduleID = CLmoduleID     [im];

    ACTS_VERBOSE(bec << " " << lydisk << " " << etamod << " " << phimod << " "
                     << side << " ");

    // Make cluster
    // TODO refactor ActsExamples::Cluster class so it is not so tedious
    Cluster cluster;

    const auto& etas = CLetas->at(im);
    const auto& phis = CLetas->at(im);
    const auto& tots = CLtots->at(im);

    const auto totalTot = std::accumulate(tots.begin(), tots.end(), 0);

    const auto [minEta, maxEta] = std::minmax_element(etas.begin(), etas.end());
    const auto [minPhi, maxPhi] = std::minmax_element(phis.begin(), phis.end());

    cluster.sizeLoc0 = *maxEta - *minEta;
    cluster.sizeLoc1 = *maxPhi - *minPhi;

    for (const auto& [eta, phi, tot] : Acts::zip(etas, phis, tots)) {
      // Make best out of what we have:
      // Weight the overall collected charge corresponding to the
      // time-over-threshold of each cell Use this as activation (does this make
      // sense?)
      auto activation = CLcharge_count[im] * tot / totalTot;

      // This bases every cluster at zero, but shouldn't matter right now
      ActsFatras::Segmentizer::Bin2D bin;
      bin[0] = eta - *minEta;
      bin[1] = phi - *minPhi;

      // Of course we have no Segment2D because this is not Fatras
      cluster.channels.emplace_back(bin, ActsFatras::Segmentizer::Segment2D{},
                                    activation);
    }

    ACTS_VERBOSE("Cluster " << im << ": " << cluster.channels.size()
                            << "cells, dimensions: " << cluster.sizeLoc0 << ", "
                            << cluster.sizeLoc1);

    clusters[im] = cluster;
  }

  // Prepare pixel space points
  SimSpacePointContainer pixelSpacePoints;

  // Prepare space-point container
  // They contain both pixel and SCT space points
  SimSpacePointContainer spacePoints;

  // Loop on space points
  for (int isp = 0; isp < nSP; isp++) {
    Acts::Vector3 globalPos{SPx[isp], SPy[isp], SPz[isp]};
    double sp_covr = SPcovr[isp];
    double sp_covz = SPcovz[isp];

    // PIX=1  STRIP = 2
    auto type = SPCL2_index[isp] == -1 ? ePixel : eStrip;

    ACTS_VERBOSE("SP:: " << type << " [" << globalPos.transpose() << "] "
                         << sp_covr << " " << sp_covz);

    boost::container::static_vector<Acts::SourceLink, 2> sLinks;

    const auto cl1Index = SPCL1_index[isp];
    assert(cl1Index >= 0 && cl1Index < nCL);

    // NOTE This of course does not produce a valid Acts-stlye geometry id, but
    // we can use it for the module map
    IndexSourceLink first(Acts::GeometryIdentifier{CLmoduleID[cl1Index]},
                          cl1Index);
    sLinks.emplace_back(first);

    if (type == eStrip) {
      const auto cl2Index = SPCL2_index[isp];
      assert(cl2Index >= 0 && cl2Index < nCL);

      // NOTE This of course does not produce a valid Acts-stlye geometry id,
      // but we can use it for the module map
      IndexSourceLink second(Acts::GeometryIdentifier{CLmoduleID[cl2Index]},
                             cl2Index);
      sLinks.emplace_back(second);
    }

    SimSpacePoint sp(globalPos, std::nullopt, sp_covr, sp_covz, std::nullopt,
                     sLinks);

    if (type == ePixel) {
      pixelSpacePoints.push_back(sp);
    }

    spacePoints.push_back(sp);
  }

  ACTS_DEBUG("Created " << pixelSpacePoints.size() << " "
                        << " pixel space points");

  ACTS_DEBUG("Created " << spacePoints.size() << " overall space points");

  m_outputPixelSpacePoints(ctx, std::move(pixelSpacePoints));
  m_outputSpacePoints(ctx, std::move(spacePoints));
  m_outputClusters(ctx, std::move(clusters));

  return ProcessCode::SUCCESS;
}

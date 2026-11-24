// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootClusterWriter.hpp"

#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <ios>

#include <TFile.h>
#include <TTree.h>

namespace ActsExamples {

RootClusterWriter::RootClusterWriter(const RootClusterWriter::Config& config,
                                     Acts::Logging::Level level)
    : WriterT(config.inputClusters, "RootClusterWriter", level), m_cfg(config) {
  // Input container for clusters is already checked by base constructor

  // Setup ROOT File
  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }

  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), "Clusters");

  m_outputTree->Branch("event_nr", &m_eventNr);
  m_outputTree->Branch("cluster_id", &m_clusterId);
  m_outputTree->Branch("volume_id", &m_volumeId);
  m_outputTree->Branch("layer_id", &m_layerId);
  m_outputTree->Branch("surface_id", &m_surfaceId);
  m_outputTree->Branch("extra_id", &m_extraId);

  m_outputTree->Branch("clus_size", &m_clusterSize);
  m_outputTree->Branch("clus_size_loc0", &m_clusterSizeLoc0);
  m_outputTree->Branch("clus_size_loc1", &m_clusterSizeLoc1);
  m_outputTree->Branch("channel_loc0", &m_channelLoc0);
  m_outputTree->Branch("channel_loc1", &m_channelLoc1);
  m_outputTree->Branch("channel_value", &m_channelValue);

  m_outputTree->Branch("clus_gx", &m_globalX);
  m_outputTree->Branch("clus_gy", &m_globalY);
  m_outputTree->Branch("clus_gz", &m_globalZ);

  // Only filled for digitization inputs that provide the cluster shape,
  // e.g. the Athena dump reader
  m_outputTree->Branch("loc_direction_x", &m_localDirectionX);
  m_outputTree->Branch("loc_direction_y", &m_localDirectionY);
  m_outputTree->Branch("loc_direction_z", &m_localDirectionZ);
  m_outputTree->Branch("length_direction_x", &m_lengthDirectionX);
  m_outputTree->Branch("length_direction_y", &m_lengthDirectionY);
  m_outputTree->Branch("length_direction_z", &m_lengthDirectionZ);
  m_outputTree->Branch("loc_eta", &m_localEta);
  m_outputTree->Branch("loc_phi", &m_localPhi);
  m_outputTree->Branch("glob_eta", &m_globalEta);
  m_outputTree->Branch("glob_phi", &m_globalPhi);
  m_outputTree->Branch("eta_angle", &m_etaAngle);
  m_outputTree->Branch("phi_angle", &m_phiAngle);
}

RootClusterWriter::~RootClusterWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ProcessCode RootClusterWriter::finalize() {
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();

  return ProcessCode::SUCCESS;
}

ProcessCode RootClusterWriter::writeT(const AlgorithmContext& ctx,
                                      const ClusterContainer& clusters) {
  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  for (Index clusterIdx = 0u; clusterIdx < clusters.size(); ++clusterIdx) {
    const Cluster& cluster = clusters[clusterIdx];

    m_eventNr = static_cast<int>(ctx.eventNumber);
    m_clusterId = clusterIdx;
    m_volumeId = static_cast<int>(cluster.geometryId.volume());
    m_layerId = static_cast<int>(cluster.geometryId.layer());
    m_surfaceId = static_cast<int>(cluster.geometryId.sensitive());
    m_extraId = static_cast<int>(cluster.geometryId.extra());

    m_clusterSize = static_cast<int>(cluster.channels.size());
    m_clusterSizeLoc0 = static_cast<int>(cluster.sizeLoc0);
    m_clusterSizeLoc1 = static_cast<int>(cluster.sizeLoc1);
    for (const auto& channel : cluster.channels) {
      m_channelLoc0.push_back(static_cast<int>(channel.bin[0]));
      m_channelLoc1.push_back(static_cast<int>(channel.bin[1]));
      m_channelValue.push_back(static_cast<float>(channel.activation));
    }

    m_globalX = static_cast<float>(cluster.globalPosition.x());
    m_globalY = static_cast<float>(cluster.globalPosition.y());
    m_globalZ = static_cast<float>(cluster.globalPosition.z());
    m_localDirectionX = static_cast<float>(cluster.localDirection.x());
    m_localDirectionY = static_cast<float>(cluster.localDirection.y());
    m_localDirectionZ = static_cast<float>(cluster.localDirection.z());
    m_lengthDirectionX = static_cast<float>(cluster.lengthDirection.x());
    m_lengthDirectionY = static_cast<float>(cluster.lengthDirection.y());
    m_lengthDirectionZ = static_cast<float>(cluster.lengthDirection.z());
    m_localEta = cluster.localEta;
    m_localPhi = cluster.localPhi;
    m_globalEta = cluster.globalEta;
    m_globalPhi = cluster.globalPhi;
    m_etaAngle = cluster.etaAngle;
    m_phiAngle = cluster.phiAngle;

    m_outputTree->Fill();

    m_channelLoc0.clear();
    m_channelLoc1.clear();
    m_channelValue.clear();
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

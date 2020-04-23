// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <mutex>

#include "ACTFW/EventData/GeometryContainers.hpp"
#include "ACTFW/Framework/WriterT.hpp"
#include "Acts/Plugins/Digitization/PlanarModuleCluster.hpp"

class TFile;
class TTree;

namespace FW {

/// @class RootPlanarClusterWriter
///
/// Write out a planar cluster collection into a root file
/// to avoid immense long vectors, each cluster is one entry
/// in the root file for optimised data writing speed
/// The event number is part of the written data.
///
/// A common file can be provided for to the writer to attach his TTree,
/// this is done by setting the Config::rootFile pointer to an existing file
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class RootPlanarClusterWriter
    : public WriterT<GeometryIdMultimap<Acts::PlanarModuleCluster>> {
 public:
  struct Config {
    /// Which cluster collection to write.
    std::string inputClusters;
    /// Which simulated (truth) hits collection to use.
    std::string inputSimulatedHits;
    std::string filePath = "";          ///< path of the output file
    std::string fileMode = "RECREATE";  ///< file access mode
    std::string treeName = "clusters";  ///< name of the output tree
    TFile* rootFile = nullptr;          ///< common root file
  };

  /// Constructor with
  /// @param cfg configuration struct
  /// @param output logging level
  RootPlanarClusterWriter(const Config& cfg, Acts::Logging::Level lvl);

  /// Virtual destructor
  ~RootPlanarClusterWriter() override;

  /// End-of-run hook
  ProcessCode endRun() final override;

 protected:
  /// This implementation holds the actual writing method
  /// and is called by the WriterT<>::write interface
  ///
  /// @param ctx The Algorithm context with per event information
  /// @param clusters is the data to be written out
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const GeometryIdMultimap<Acts::PlanarModuleCluster>&
                         clusters) final override;

 private:
  Config m_cfg;                    ///< the configuration object
  std::mutex m_writeMutex;         ///< protect multi-threaded writes
  TFile* m_outputFile;             ///< the output file
  TTree* m_outputTree;             ///< the output tree
  int m_eventNr;                   ///< the event number of
  int m_volumeID;                  ///< volume identifier
  int m_layerID;                   ///< layer identifier
  int m_surfaceID;                 ///< surface identifier
  float m_x;                       ///< global x
  float m_y;                       ///< global y
  float m_z;                       ///< global z
  float m_t;                       ///< global t
  float m_lx;                      ///< local lx
  float m_ly;                      ///< local ly
  float m_cov_lx;                  ///< local covariance lx
  float m_cov_ly;                  ///< local covariance ly
  std::vector<int> m_cell_IDx;     ///< cell ID in lx
  std::vector<int> m_cell_IDy;     ///< cell ID in ly
  std::vector<float> m_cell_lx;    ///< local cell position x
  std::vector<float> m_cell_ly;    ///< local cell position y
  std::vector<float> m_cell_data;  ///< local cell position y

  // (optional) the truth position
  std::vector<float> m_t_gx;  ///< truth position global x
  std::vector<float> m_t_gy;  ///< truth position global y
  std::vector<float> m_t_gz;  ///< truth position global z
  std::vector<float> m_t_gt;  ///< truth time t
  std::vector<float> m_t_lx;  ///< truth position local x
  std::vector<float> m_t_ly;  ///< truth position local y
  std::vector<unsigned long>
      m_t_barcode;  ///< associated truth particle barcode
};

}  // namespace FW

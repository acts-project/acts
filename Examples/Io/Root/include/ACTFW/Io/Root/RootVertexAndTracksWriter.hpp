// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <mutex>

#include "ACTFW/Framework/WriterT.hpp"
#include "ACTFW/TruthTracking/VertexAndTracks.hpp"

class TFile;
class TTree;

namespace FW {

/// Write out vertices together with associated tracks into a TTree
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
///
/// A common file can be provided for to the writer to attach his TTree,
/// this is done by setting the Config::rootFile pointer to an existing file
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class RootVertexAndTracksWriter final
    : public WriterT<std::vector<VertexAndTracks>> {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::string collection;             ///< particle collection to write
    std::string filePath;               ///< path of the output file
    std::string fileMode = "RECREATE";  ///< file access mode
    std::string treeName = "event";     ///< name of the output tree
    TFile* rootFile = nullptr;          ///< common root file
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  RootVertexAndTracksWriter(const Config& cfg, Acts::Logging::Level lvl);

  /// Virtual destructor
  ~RootVertexAndTracksWriter() final override;

  /// End-of-run hook
  ProcessCode endRun() final override;

 protected:
  /// @brief Write method called by the base class
  /// @param [in] context is the algorithm context for event information
  /// @param [in] vertexAndTracksCollection is the VertexAndTracks collection
  ProcessCode writeT(const AlgorithmContext& context,
                     const std::vector<VertexAndTracks>&
                         vertexAndTracksCollection) final override;

 private:
  Config m_cfg;             ///< The config class
  std::mutex m_writeMutex;  ///< Mutex used to protect multi-threaded writes
  TFile* m_outputFile{nullptr};  ///< The output file
  TTree* m_outputTree{nullptr};  ///< The output tree
  int m_eventNr{0};              ///< the event number of

  /// The vertex positions
  std::vector<double> m_vx;
  std::vector<double> m_vy;
  std::vector<double> m_vz;

  /// The track parameter
  std::vector<double> m_d0;
  std::vector<double> m_z0;
  std::vector<double> m_phi;
  std::vector<double> m_theta;
  std::vector<double> m_qp;
  std::vector<double> m_time;
  std::vector<int> m_vtxID;

  /// The track covariance matrix
  std::vector<double> m_cov11;
  std::vector<double> m_cov12;
  std::vector<double> m_cov13;
  std::vector<double> m_cov14;
  std::vector<double> m_cov15;
  std::vector<double> m_cov16;

  std::vector<double> m_cov21;
  std::vector<double> m_cov22;
  std::vector<double> m_cov23;
  std::vector<double> m_cov24;
  std::vector<double> m_cov25;
  std::vector<double> m_cov26;

  std::vector<double> m_cov31;
  std::vector<double> m_cov32;
  std::vector<double> m_cov33;
  std::vector<double> m_cov34;
  std::vector<double> m_cov35;
  std::vector<double> m_cov36;

  std::vector<double> m_cov41;
  std::vector<double> m_cov42;
  std::vector<double> m_cov43;
  std::vector<double> m_cov44;
  std::vector<double> m_cov45;
  std::vector<double> m_cov46;

  std::vector<double> m_cov51;
  std::vector<double> m_cov52;
  std::vector<double> m_cov53;
  std::vector<double> m_cov54;
  std::vector<double> m_cov55;
  std::vector<double> m_cov56;

  std::vector<double> m_cov61;
  std::vector<double> m_cov62;
  std::vector<double> m_cov63;
  std::vector<double> m_cov64;
  std::vector<double> m_cov65;
  std::vector<double> m_cov66;

  /// Pointers to the vectors
  std::vector<double>* m_ptrVx = &m_vx;
  std::vector<double>* m_ptrVy = &m_vy;
  std::vector<double>* m_ptrVz = &m_vz;
  std::vector<double>* m_ptrD0 = &m_d0;
  std::vector<double>* m_ptrZ0 = &m_z0;
  std::vector<double>* m_ptrPhi = &m_phi;
  std::vector<double>* m_ptrTheta = &m_theta;
  std::vector<double>* m_ptrQP = &m_qp;
  std::vector<double>* m_ptrTime = &m_time;
  std::vector<int>* m_ptrVtxID = &m_vtxID;

  std::vector<double>* m_ptrCov11 = &m_cov11;
  std::vector<double>* m_ptrCov12 = &m_cov12;
  std::vector<double>* m_ptrCov13 = &m_cov13;
  std::vector<double>* m_ptrCov14 = &m_cov14;
  std::vector<double>* m_ptrCov15 = &m_cov15;
  std::vector<double>* m_ptrCov16 = &m_cov16;

  std::vector<double>* m_ptrCov21 = &m_cov21;
  std::vector<double>* m_ptrCov22 = &m_cov22;
  std::vector<double>* m_ptrCov23 = &m_cov23;
  std::vector<double>* m_ptrCov24 = &m_cov24;
  std::vector<double>* m_ptrCov25 = &m_cov25;
  std::vector<double>* m_ptrCov26 = &m_cov26;

  std::vector<double>* m_ptrCov31 = &m_cov31;
  std::vector<double>* m_ptrCov32 = &m_cov32;
  std::vector<double>* m_ptrCov33 = &m_cov33;
  std::vector<double>* m_ptrCov34 = &m_cov34;
  std::vector<double>* m_ptrCov35 = &m_cov35;
  std::vector<double>* m_ptrCov36 = &m_cov36;

  std::vector<double>* m_ptrCov41 = &m_cov41;
  std::vector<double>* m_ptrCov42 = &m_cov42;
  std::vector<double>* m_ptrCov43 = &m_cov43;
  std::vector<double>* m_ptrCov44 = &m_cov44;
  std::vector<double>* m_ptrCov45 = &m_cov45;
  std::vector<double>* m_ptrCov46 = &m_cov46;

  std::vector<double>* m_ptrCov51 = &m_cov51;
  std::vector<double>* m_ptrCov52 = &m_cov52;
  std::vector<double>* m_ptrCov53 = &m_cov53;
  std::vector<double>* m_ptrCov54 = &m_cov54;
  std::vector<double>* m_ptrCov55 = &m_cov55;
  std::vector<double>* m_ptrCov56 = &m_cov56;

  std::vector<double>* m_ptrCov61 = &m_cov61;
  std::vector<double>* m_ptrCov62 = &m_cov62;
  std::vector<double>* m_ptrCov63 = &m_cov63;
  std::vector<double>* m_ptrCov64 = &m_cov64;
  std::vector<double>* m_ptrCov65 = &m_cov65;
  std::vector<double>* m_ptrCov66 = &m_cov66;

  /// @brief Clears all vectors
  void ClearAll();
};

}  // namespace FW

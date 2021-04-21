// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <mutex>
#include <type_traits>

#include "ACTFW/Framework/WriterT.hpp"
#include "ACTFW/TruthTracking/VertexAndTracks.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingError.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"

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
class RootRecVertexWriter final
    : public WriterT<std::vector<Acts::Vertex<Acts::BoundParameters>>> {
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
  RootRecVertexWriter(const Config& cfg, Acts::Logging::Level lvl);

  /// Virtual destructor
  ~RootRecVertexWriter() final override;

  /// End-of-run hook
  ProcessCode endRun() final override;

 protected:
  /// @brief Write method called by the base class
  /// @param [in] context is the algorithm context for event information
  /// @param [in] vertexCollection is the RecVertex collection // Haoran
  /// modified. See AdaptiveMultiVertexFinderAlgorithm.cpp
  ProcessCode writeT(const AlgorithmContext& context,
                     const std::vector<Acts::Vertex<Acts::BoundParameters>>&
                         vertexCollection) final override;

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

  // std::vector<std::pair<double, double>> m_vtx_fitquality;
  std::vector<double> m_vtx_fitquality_chiSquared;
  std::vector<double> m_vtx_fitquality_numberDoF;

  /// The vertex covariance matrix
  std::vector<double> m_vtx_cov11;
  std::vector<double> m_vtx_cov12;
  std::vector<double> m_vtx_cov13;
  std::vector<double> m_vtx_cov14;

  std::vector<double> m_vtx_cov21;
  std::vector<double> m_vtx_cov22;
  std::vector<double> m_vtx_cov23;
  std::vector<double> m_vtx_cov24;

  std::vector<double> m_vtx_cov31;
  std::vector<double> m_vtx_cov32;
  std::vector<double> m_vtx_cov33;
  std::vector<double> m_vtx_cov34;

  std::vector<double> m_vtx_cov41;
  std::vector<double> m_vtx_cov42;
  std::vector<double> m_vtx_cov43;
  std::vector<double> m_vtx_cov44;

  /// The track parameter
  std::vector<double> m_d0;
  std::vector<double> m_z0;
  std::vector<double> m_phi;
  std::vector<double> m_theta;
  std::vector<double> m_qp;
  std::vector<double> m_time;
  std::vector<int> m_vtxID;

  /// The track covariance matrix
  std::vector<double> m_trk_cov11;
  std::vector<double> m_trk_cov12;
  std::vector<double> m_trk_cov13;
  std::vector<double> m_trk_cov14;
  std::vector<double> m_trk_cov15;
  std::vector<double> m_trk_cov16;

  std::vector<double> m_trk_cov21;
  std::vector<double> m_trk_cov22;
  std::vector<double> m_trk_cov23;
  std::vector<double> m_trk_cov24;
  std::vector<double> m_trk_cov25;
  std::vector<double> m_trk_cov26;

  std::vector<double> m_trk_cov31;
  std::vector<double> m_trk_cov32;
  std::vector<double> m_trk_cov33;
  std::vector<double> m_trk_cov34;
  std::vector<double> m_trk_cov35;
  std::vector<double> m_trk_cov36;

  std::vector<double> m_trk_cov41;
  std::vector<double> m_trk_cov42;
  std::vector<double> m_trk_cov43;
  std::vector<double> m_trk_cov44;
  std::vector<double> m_trk_cov45;
  std::vector<double> m_trk_cov46;

  std::vector<double> m_trk_cov51;
  std::vector<double> m_trk_cov52;
  std::vector<double> m_trk_cov53;
  std::vector<double> m_trk_cov54;
  std::vector<double> m_trk_cov55;
  std::vector<double> m_trk_cov56;

  std::vector<double> m_trk_cov61;
  std::vector<double> m_trk_cov62;
  std::vector<double> m_trk_cov63;
  std::vector<double> m_trk_cov64;
  std::vector<double> m_trk_cov65;
  std::vector<double> m_trk_cov66;

  /// Pointers to the vectors
  std::vector<double>* m_ptrVx = &m_vx;
  std::vector<double>* m_ptrVy = &m_vy;
  std::vector<double>* m_ptrVz = &m_vz;
  // std::vector<std::pair<double, double>>* m_ptr_vtx_fitquality =
  //     &m_vtx_fitquality;

  std::vector<double>* m_ptr_vtx_fitquality_chiSquared =
      &m_vtx_fitquality_chiSquared;
  std::vector<double>* m_ptr_vtx_fitquality_numberDoF =
      &m_vtx_fitquality_numberDoF;

  std::vector<double>* m_ptrD0 = &m_d0;
  std::vector<double>* m_ptrZ0 = &m_z0;
  std::vector<double>* m_ptrPhi = &m_phi;
  std::vector<double>* m_ptrTheta = &m_theta;
  std::vector<double>* m_ptrQP = &m_qp;
  std::vector<double>* m_ptrTime = &m_time;
  std::vector<int>* m_ptrVtxID = &m_vtxID;

  std::vector<double>* m_ptr_vtx_Cov11 = &m_vtx_cov11;
  std::vector<double>* m_ptr_vtx_Cov12 = &m_vtx_cov12;
  std::vector<double>* m_ptr_vtx_Cov13 = &m_vtx_cov13;
  std::vector<double>* m_ptr_vtx_Cov14 = &m_vtx_cov14;

  std::vector<double>* m_ptr_vtx_Cov21 = &m_vtx_cov21;
  std::vector<double>* m_ptr_vtx_Cov22 = &m_vtx_cov22;
  std::vector<double>* m_ptr_vtx_Cov23 = &m_vtx_cov23;
  std::vector<double>* m_ptr_vtx_Cov24 = &m_vtx_cov24;

  std::vector<double>* m_ptr_vtx_Cov31 = &m_vtx_cov31;
  std::vector<double>* m_ptr_vtx_Cov32 = &m_vtx_cov32;
  std::vector<double>* m_ptr_vtx_Cov33 = &m_vtx_cov33;
  std::vector<double>* m_ptr_vtx_Cov34 = &m_vtx_cov34;

  std::vector<double>* m_ptr_vtx_Cov41 = &m_vtx_cov41;
  std::vector<double>* m_ptr_vtx_Cov42 = &m_vtx_cov42;
  std::vector<double>* m_ptr_vtx_Cov43 = &m_vtx_cov43;
  std::vector<double>* m_ptr_vtx_Cov44 = &m_vtx_cov44;

  std::vector<double>* m_ptr_trk_Cov11 = &m_trk_cov11;
  std::vector<double>* m_ptr_trk_Cov12 = &m_trk_cov12;
  std::vector<double>* m_ptr_trk_Cov13 = &m_trk_cov13;
  std::vector<double>* m_ptr_trk_Cov14 = &m_trk_cov14;
  std::vector<double>* m_ptr_trk_Cov15 = &m_trk_cov15;
  std::vector<double>* m_ptr_trk_Cov16 = &m_trk_cov16;

  std::vector<double>* m_ptr_trk_Cov21 = &m_trk_cov21;
  std::vector<double>* m_ptr_trk_Cov22 = &m_trk_cov22;
  std::vector<double>* m_ptr_trk_Cov23 = &m_trk_cov23;
  std::vector<double>* m_ptr_trk_Cov24 = &m_trk_cov24;
  std::vector<double>* m_ptr_trk_Cov25 = &m_trk_cov25;
  std::vector<double>* m_ptr_trk_Cov26 = &m_trk_cov26;

  std::vector<double>* m_ptr_trk_Cov31 = &m_trk_cov31;
  std::vector<double>* m_ptr_trk_Cov32 = &m_trk_cov32;
  std::vector<double>* m_ptr_trk_Cov33 = &m_trk_cov33;
  std::vector<double>* m_ptr_trk_Cov34 = &m_trk_cov34;
  std::vector<double>* m_ptr_trk_Cov35 = &m_trk_cov35;
  std::vector<double>* m_ptr_trk_Cov36 = &m_trk_cov36;

  std::vector<double>* m_ptr_trk_Cov41 = &m_trk_cov41;
  std::vector<double>* m_ptr_trk_Cov42 = &m_trk_cov42;
  std::vector<double>* m_ptr_trk_Cov43 = &m_trk_cov43;
  std::vector<double>* m_ptr_trk_Cov44 = &m_trk_cov44;
  std::vector<double>* m_ptr_trk_Cov45 = &m_trk_cov45;
  std::vector<double>* m_ptr_trk_Cov46 = &m_trk_cov46;

  std::vector<double>* m_ptr_trk_Cov51 = &m_trk_cov51;
  std::vector<double>* m_ptr_trk_Cov52 = &m_trk_cov52;
  std::vector<double>* m_ptr_trk_Cov53 = &m_trk_cov53;
  std::vector<double>* m_ptr_trk_Cov54 = &m_trk_cov54;
  std::vector<double>* m_ptr_trk_Cov55 = &m_trk_cov55;
  std::vector<double>* m_ptr_trk_Cov56 = &m_trk_cov56;

  std::vector<double>* m_ptr_trk_Cov61 = &m_trk_cov61;
  std::vector<double>* m_ptr_trk_Cov62 = &m_trk_cov62;
  std::vector<double>* m_ptr_trk_Cov63 = &m_trk_cov63;
  std::vector<double>* m_ptr_trk_Cov64 = &m_trk_cov64;
  std::vector<double>* m_ptr_trk_Cov65 = &m_trk_cov65;
  std::vector<double>* m_ptr_trk_Cov66 = &m_trk_cov66;

  /// @brief Clears all vectors
  void ClearAll();
};

}  // namespace FW

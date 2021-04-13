// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "Acts/Vertexing/Vertex.hpp"

#include <mutex>
#include <vector>

class TFile;
class TTree;

namespace ActsExamples {

/// @class RootVertexPerformanceWriter
class RootVertexPerformanceWriter final
    : public WriterT<std::vector<Acts::Vertex<Acts::BoundTrackParameters>>> {
 public:
  struct Config {
    /// Input truth trajectories collection
    std::string inputTruthParticles;
    /// Input vertex collection.
    std::string inputVertices;
    /// output directory.
    std::string outputDir;
    /// output filename.
    std::string outputFilename = "vertexingperformance.root";
    /// name of the output tree.
    std::string outputTreename = "vertextree";
    /// file access mode.
    std::string fileMode = "RECREATE";
    /// common root file.
    TFile* rootFile = nullptr;
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param level Message level declaration
  RootVertexPerformanceWriter(const Config& cfg, Acts::Logging::Level lvl);
  ~RootVertexPerformanceWriter() final override;

  /// End-of-run hook
  ProcessCode endRun() final override;

 protected:
  /// @brief Write method called by the base class
  /// @param [in] ctx is the algorithm context for event information
  ProcessCode writeT(const AlgorithmContext& ctx, const std::vector<Acts::Vertex<Acts::BoundTrackParameters>>& vertices) final override;

 private:
  Config m_cfg;             ///< The config class
  std::mutex m_writeMutex;  ///< Mutex used to protect multi-threaded writes
  TFile* m_outputFile{nullptr};   ///< The output file
  TTree* m_outputTree{nullptr};   ///< The output tree

  std::vector<float> m_diffx;
  std::vector<float> m_diffy;
  std::vector<float> m_diffz;

  int m_nrecoVtx = -1;
  int m_ntrueVtx = -1;
};

}  // namespace ActsExamples

// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <mutex>
#include <vector>

class TFile;
class TTree;

namespace ActsExamples {

/// @class RootVertexPerformanceWriter
///
/// Writes out the number of reconstructed primary vertices along with
/// the number of primary vertices in detector acceptance as well as
/// reconstructable primary vertices after track fitting.
/// Additionally it matches the reco vertices to their truth vertices
/// and write out the difference in x,y and z position.
class RootVertexPerformanceWriter final
    : public WriterT<std::vector<Acts::Vertex<Acts::BoundTrackParameters>>> {
 public:
  struct Config {
    /// All input truth particle collection
    std::string inputAllTruthParticles;
    /// Selected input truth particle collection
    std::string inputSelectedTruthParticles;
    /// Truth particles associated to fitted tracks
    std::string inputAssociatedTruthParticles;
    /// Selected fitted tracks
    std::string inputFittedTracks;
    /// Selected fitted tracks indices (points from `inputSelectedFittedTracks`
    /// to `inputAllFittedTracks`). If empty all we assume selected == all.
    std::string inputFittedTracksIndices;
    /// All event fitted tracks tips (points from `inputAllFittedTracks` to
    /// `inputTrajectories`)
    std::string inputAllFittedTracksTips;
    /// Trajectories object from track finidng
    std::string inputTrajectories;
    /// Input hit-particles map collection.
    std::string inputMeasurementParticlesMap;
    /// Input vertex collection.
    std::string inputVertices;
    /// Input reconstruction time.
    std::string inputTime;
    /// Output filename.
    std::string filePath = "vertexingperformance.root";
    /// Name of the output tree.
    std::string treeName = "vertextree";
    /// File access mode.
    std::string fileMode = "RECREATE";
    /// Minimum fraction of tracks matched between truth
    /// and reco vertices to be matched for resolution plots
    double minTrackVtxMatchFraction = 0.5;
    /// Minimum fraction of hits associated to particle to consider
    /// as truth matched
    double truthMatchProbMin = 0.5;
  };

  /// Constructor
  ///
  /// @param config Configuration struct
  /// @param level Message level declaration
  RootVertexPerformanceWriter(const Config& config, Acts::Logging::Level level);

  ~RootVertexPerformanceWriter() final override;

  /// End-of-run hook
  ProcessCode endRun() final override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 protected:
  /// @brief Write method called by the base class
  /// @param [in] ctx is the algorithm context for event information
  ProcessCode writeT(
      const AlgorithmContext& ctx,
      const std::vector<Acts::Vertex<Acts::BoundTrackParameters>>& vertices)
      final override;

 private:
  Config m_cfg;             ///< The config class
  std::mutex m_writeMutex;  ///< Mutex used to protect multi-threaded writes
  TFile* m_outputFile{nullptr};  ///< The output file
  TTree* m_outputTree{nullptr};  ///< The output tree

  std::vector<float>
      m_diffx;  ///< Difference in x positon between reco and true vtx
  std::vector<float>
      m_diffy;  ///< Difference in y positon between reco and true vtx
  std::vector<float>
      m_diffz;  ///< Difference in z positon between reco and true vtx

  std::vector<float> m_truthX;
  std::vector<float> m_truthY;
  std::vector<float> m_truthZ;

  std::vector<float> m_recoX;
  std::vector<float> m_recoY;
  std::vector<float> m_recoZ;

  std::vector<float> m_covXX;
  std::vector<float> m_covYY;
  std::vector<float> m_covXY;
  std::vector<float> m_covYX;
  std::vector<float> m_trackVtxMatchFraction;

  int m_nrecoVtx = -1;           ///< Number of reconstructed vertices
  int m_ntrueVtx = -1;           ///< Number of true vertices
  int m_nVtxDetAcceptance = -1;  ///< Number of vertices in detector acceptance
  int m_nVtxReconstructable =
      -1;  ///< Max. number of reconstructable vertices (detector acceptance +
           ///< tracking efficiency)
  int m_timeMS = -1;  ///< Reconstruction time in ms

  int getNumberOfReconstructableVertices(
      const SimParticleContainer& collection) const;

  int getNumberOfTruePriVertices(const SimParticleContainer& collection) const;
};

}  // namespace ActsExamples

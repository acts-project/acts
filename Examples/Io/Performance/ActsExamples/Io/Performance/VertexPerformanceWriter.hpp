// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <mutex>
#include <string>
#include <vector>

class TFile;
class TTree;
namespace ActsFatras {
class Barcode;
}  // namespace ActsFatras

namespace ActsExamples {
struct AlgorithmContext;

/// @class VertexPerformanceWriter
///
/// Writes out the number of reconstructed primary vertices along with
/// the number of primary vertices in detector acceptance as well as
/// reconstructable primary vertices after track fitting.
/// Additionally it matches the reco vertices to their truth vertices
/// and write out the difference in x,y and z position.
class VertexPerformanceWriter final
    : public WriterT<std::vector<Acts::Vertex<Acts::BoundTrackParameters>>> {
 public:
  using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;

  struct Config {
    /// All input truth particle collection.
    std::string inputAllTruthParticles;
    /// Selected input truth particle collection.
    std::string inputSelectedTruthParticles;
    /// Tracks object from track finidng.
    std::string inputTracks;
    /// Optional. Truth particles associated to tracks. Using 1:1 matching if
    /// given.
    std::string inputAssociatedTruthParticles;
    /// Input hit-particles map collection.
    std::string inputMeasurementParticlesMap;
    /// Input vertex collection.
    std::string inputVertices;
    /// Magnetic field
    std::shared_ptr<Acts::MagneticFieldProvider> bField;
    /// Output filename.
    std::string filePath = "vertexingperformance.root";
    /// Name of the output tree.
    std::string treeName = "vertextree";
    /// File access mode.
    std::string fileMode = "RECREATE";
    /// Minimum fraction of tracks matched between truth
    /// and reco vertices to be matched for resolution plots.
    double minTrackVtxMatchFraction = 0.5;
    /// Minimum fraction of hits associated to particle to consider
    /// as truth matched.
    double truthMatchProbMin = 0.5;
    /// Whether information about tracks is available
    bool useTracks = true;
  };

  /// Constructor
  ///
  /// @param config Configuration struct
  /// @param level Message level declaration
  VertexPerformanceWriter(const Config& config, Acts::Logging::Level level);

  ~VertexPerformanceWriter() override;

  /// End-of-run hook
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 protected:
  /// @brief Write method called by the base class
  /// @param [in] ctx is the algorithm context for event information
  ProcessCode writeT(
      const AlgorithmContext& ctx,
      const std::vector<Acts::Vertex<Acts::BoundTrackParameters>>& vertices)
      override;

 private:
  Config m_cfg;             ///< The config class
  std::mutex m_writeMutex;  ///< Mutex used to protect multi-threaded writes
  TFile* m_outputFile{nullptr};  ///< The output file
  TTree* m_outputTree{nullptr};  ///< The output tree

  // True 4D vertex position
  std::vector<double> m_truthX;
  std::vector<double> m_truthY;
  std::vector<double> m_truthZ;
  std::vector<double> m_truthT;

  // Reconstructed 4D vertex position
  std::vector<double> m_recoX;
  std::vector<double> m_recoY;
  std::vector<double> m_recoZ;
  std::vector<double> m_recoT;

  /// Difference of reconstructed and true vertex 4D position
  std::vector<double> m_resX;
  std::vector<double> m_resY;
  std::vector<double> m_resZ;
  std::vector<double> m_resT;

  // pull(X) = (X_reco - X_true)/Var(X_reco)^(1/2)
  std::vector<double> m_pullX;
  std::vector<double> m_pullY;
  std::vector<double> m_pullZ;
  std::vector<double> m_pullT;

  // Vertex covariance
  std::vector<double> m_covXX;
  std::vector<double> m_covYY;
  std::vector<double> m_covZZ;
  std::vector<double> m_covTT;
  std::vector<double> m_covXY;
  std::vector<double> m_covXZ;
  std::vector<double> m_covXT;
  std::vector<double> m_covYZ;
  std::vector<double> m_covYT;
  std::vector<double> m_covZT;

  //--------------------------------------------------------------
  // Track-related variables are contained in a vector of vectors: The inner
  // vectors contain the values of all tracks corresponding to one vertex. The
  // outer vector can then have the same length as the flat vectors of
  // vertex-related variables (see above). E.g.,
  // m_truthPhi = ((truthPhi of 1st trk belonging to vtx 1,
  //                truthPhi of 2nd trk belonging to vtx 1, ...),
  //               (truthPhi of 1st trk belonging to vtx 2,
  //                truthPhi of 2nd trk belonging to vtx 2, ...),
  //                ...)
  //
  // True track momenta at the vertex
  std::vector<std::vector<double>> m_truthPhi;
  std::vector<std::vector<double>> m_truthTheta;
  std::vector<std::vector<double>> m_truthQOverP;

  // Reconstructed track momenta at the vertex before and after the vertex fit
  std::vector<std::vector<double>> m_recoPhi;
  std::vector<std::vector<double>> m_recoPhiFitted;
  std::vector<std::vector<double>> m_recoTheta;
  std::vector<std::vector<double>> m_recoThetaFitted;
  std::vector<std::vector<double>> m_recoQOverP;
  std::vector<std::vector<double>> m_recoQOverPFitted;

  // Difference between reconstructed momenta and true momenta
  std::vector<std::vector<double>> m_resPhi;
  std::vector<std::vector<double>> m_resPhiFitted;
  std::vector<std::vector<double>> m_resTheta;
  std::vector<std::vector<double>> m_resThetaFitted;
  std::vector<std::vector<double>> m_resQOverP;
  std::vector<std::vector<double>> m_resQOverPFitted;
  std::vector<std::vector<double>> m_momOverlap;
  std::vector<std::vector<double>> m_momOverlapFitted;

  // Pulls
  std::vector<std::vector<double>> m_pullPhi;
  std::vector<std::vector<double>> m_pullPhiFitted;
  std::vector<std::vector<double>> m_pullTheta;
  std::vector<std::vector<double>> m_pullThetaFitted;
  std::vector<std::vector<double>> m_pullQOverP;
  std::vector<std::vector<double>> m_pullQOverPFitted;

  // Number of tracks associated with truth/reconstructed vertex
  std::vector<int> m_nTracksOnTruthVertex;
  std::vector<int> m_nTracksOnRecoVertex;

  std::vector<double> m_trackVtxMatchFraction;

  /// Number of reconstructed vertices
  int m_nRecoVtx = -1;
  /// Number of true vertices
  int m_nTrueVtx = -1;
  /// Number of vertices in detector acceptance
  int m_nVtxDetAcceptance = -1;
  /// Max. number of reconstructable vertices (detector acceptance + tracking
  /// efficiency)
  int m_nVtxReconstructable = -1;

  int getNumberOfReconstructableVertices(
      const SimParticleContainer& collection) const;

  int getNumberOfTruePriVertices(const SimParticleContainer& collection) const;

  ReadDataHandle<SimParticleContainer> m_inputAllTruthParticles{
      this, "InputAllTruthParticles"};

  ReadDataHandle<SimParticleContainer> m_inputSelectedTruthParticles{
      this, "InputSelectedTruthParticles"};

  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "InputTracks"};

  ReadDataHandle<SimParticleContainer> m_inputAssociatedTruthParticles{
      this, "InputAssociatedTruthParticles"};

  ReadDataHandle<HitParticlesMap> m_inputMeasurementParticlesMap{
      this, "InputMeasurementParticlesMap"};
};

}  // namespace ActsExamples

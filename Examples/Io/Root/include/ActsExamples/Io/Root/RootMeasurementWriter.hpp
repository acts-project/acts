// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsFatras/Digitization/Channelizer.hpp"

#include <array>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

#include <TTree.h>

class TFile;
class TTree;
namespace Acts {
class Surface;
class TrackingGeometry;
}  // namespace Acts

namespace ActsExamples {
struct AlgorithmContext;

/// @class RootMeasurementWriter
///
/// Write out a planar cluster collection into a root file
/// to avoid immense long vectors, each cluster is one entry
/// in the root file for optimised data writing speed
/// The event number is part of the written data.
///
/// A common file can be provided for the writer to attach his TTree,
/// this is done by setting the Config::rootFile pointer to an existing file
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class RootMeasurementWriter final : public WriterT<MeasurementContainer> {
 public:
  struct Config {
    /// Which measurement collection to write.
    std::string inputMeasurements;
    /// Which cluster collection to write (optional)
    std::string inputClusters;
    /// Which simulated (truth) hits collection to use.
    std::string inputSimHits;
    /// Input collection to map measured hits to simulated hits.
    std::string inputMeasurementSimHitsMap;
    std::string filePath = "";          ///< path of the output file
    std::string fileMode = "RECREATE";  ///< file access mode
    /// The indices for this digitization configurations
    Acts::GeometryHierarchyMap<std::vector<Acts::BoundIndices>> boundIndices;
    /// Map of the geometry identifier to the surface
    std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*>
        surfaceByIdentifier;
  };

  struct DigitizationTree {
    const std::array<std::string, Acts::eBoundSize> bNames = {
        "loc0", "loc1", "phi", "theta", "qop", "time"};

    TTree* tree = nullptr;
    // Identification parameters
    int eventNr = 0;
    int volumeID = 0;
    int layerID = 0;
    int surfaceID = 0;

    /// Type 0 - free, 1 - bound
    int measType = 1;

    /// Truth parameters
    float trueBound[Acts::eBoundSize] = {};
    float trueGx = 0.;
    float trueGy = 0.;
    float trueGz = 0.;
    float incidentPhi = 0.;
    float incidentTheta = 0.;

    /// Reconstruction information
    float recBound[Acts::eBoundSize] = {};
    float varBound[Acts::eBoundSize] = {};

    /// Cluster information comprised of
    /// nch :  number of channels
    /// cSize : cluster size in loc0 and loc1
    /// chId : channel identification
    /// chValue: value/activation of the channel
    int nch = 0;
    int cSize[2] = {};
    std::array<std::vector<int>*, 2> chId = {nullptr, nullptr};
    std::vector<float>* chValue = nullptr;

    /// Setup helper to create the tree and
    /// register the branches
    ///
    /// @param treeName the name of the tree to be registered
    void setupTree(const std::string& treeName) {
      tree = new TTree(treeName.c_str(), treeName.c_str());
      // Declare the branches
      tree->Branch("event_nr", &eventNr);
      tree->Branch("volume_id", &volumeID);
      tree->Branch("layer_id", &layerID);
      tree->Branch("surface_id", &surfaceID);
      tree->Branch("measurement_type", &measType);
      for (unsigned int ib = 0; ib < Acts::eBoundSize; ++ib) {
        if (ib != Acts::eBoundQOverP) {
          tree->Branch(("true_" + bNames[ib]).c_str(), &trueBound[ib]);
        }
      }
      tree->Branch("true_x", &trueGx);
      tree->Branch("true_y", &trueGy);
      tree->Branch("true_z", &trueGz);
      tree->Branch("true_incident_phi", &incidentPhi);
      tree->Branch("true_incident_theta", &incidentTheta);
    }

    /// Constructor from GeometryIdentifier
    DigitizationTree(Acts::GeometryIdentifier geoID) {
      auto vID = geoID.volume();
      auto lID = geoID.layer();
      auto mID = geoID.sensitive();
      std::string treeName = "vol" + std::to_string(vID);
      if (lID > 0) {
        treeName += "_lay" + std::to_string(lID);
      }
      if (mID > 0) {
        treeName += "_mod" + std::to_string(mID);
      }
      setupTree(treeName);
    }

    /// Non-trivial destructor for memory cleanup
    ~DigitizationTree() {
      delete chId[0];
      delete chId[1];
      delete chValue;
    }

    /// Setup the dimension depended branches
    ///
    /// @param i the bound index in question
    void setupBoundRecBranch(Acts::BoundIndices i) {
      tree->Branch(("rec_" + bNames[i]).c_str(), &recBound[i]);
      tree->Branch(("var_" + bNames[i]).c_str(), &varBound[i]);
    }

    /// Setup the cluster related branch
    ///
    /// @param bIndices the bound indices to be written
    void setupClusterBranch(const std::vector<Acts::BoundIndices>& bIndices) {
      chValue = new std::vector<float>;
      tree->Branch("clus_size", &nch);
      tree->Branch("channel_value", &chValue);
      // Both are allocated, but only relevant ones are set
      chId[0] = new std::vector<int>;
      chId[1] = new std::vector<int>;
      for (const auto& ib : bIndices) {
        if (static_cast<unsigned int>(ib) < 2) {
          tree->Branch(("channel_" + bNames[ib]).c_str(), &chId[ib]);
          tree->Branch(("clus_size_" + bNames[ib]).c_str(), &cSize[ib]);
        }
      }
    }

    /// Convenience function to register idenfication
    ///
    /// @param eventNr The event number
    /// @param geoID The geometry identifier of the measurement
    void fillIdentification(int evnt, Acts::GeometryIdentifier geoId) {
      eventNr = evnt;
      volumeID = geoId.volume();
      layerID = geoId.layer();
      surfaceID = geoId.sensitive();
    }

    /// Convenience function to register the truth parameters
    ///
    /// @param lp The true local position
    /// @param xt The true 4D global position
    /// @param dir The true particle direction
    void fillTruthParameters(const Acts::Vector2& lp, const Acts::Vector4& xt,
                             const Acts::Vector3& dir,
                             const std::pair<double, double> angles) {
      trueBound[Acts::eBoundLoc0] = lp[Acts::eBoundLoc0];
      trueBound[Acts::eBoundLoc1] = lp[Acts::eBoundLoc1];
      trueBound[Acts::eBoundPhi] = Acts::VectorHelpers::phi(dir);
      trueBound[Acts::eBoundTheta] = Acts::VectorHelpers::theta(dir);
      trueBound[Acts::eBoundTime] = xt[Acts::eTime];

      trueGx = xt[Acts::ePos0];
      trueGy = xt[Acts::ePos1];
      trueGz = xt[Acts::ePos2];

      incidentPhi = angles.first;
      incidentTheta = angles.second;
    }

    /// Convenience function to fill bound parameters
    ///
    /// @tparam measurement_t Type of the parameter set
    ///
    /// @param m The measurement set
    template <typename measurement_t>
    void fillBoundMeasurement(const measurement_t& m) {
      Acts::BoundVector fullVect = m.expander() * m.parameters();
      recBound[Acts::eBoundLoc0] = fullVect[Acts::eBoundLoc0];
      recBound[Acts::eBoundLoc1] = fullVect[Acts::eBoundLoc1];
      recBound[Acts::eBoundPhi] = fullVect[Acts::eBoundPhi];
      recBound[Acts::eBoundTheta] = fullVect[Acts::eBoundTheta];
      recBound[Acts::eBoundTime] = fullVect[Acts::eBoundTime];

      Acts::BoundSquareMatrix fullVar =
          m.expander() * m.covariance() * m.expander().transpose();
      varBound[Acts::eBoundLoc0] = fullVar(Acts::eBoundLoc0, Acts::eBoundLoc0);
      varBound[Acts::eBoundLoc1] = fullVar(Acts::eBoundLoc1, Acts::eBoundLoc1);
      varBound[Acts::eBoundPhi] = fullVar(Acts::eBoundPhi, Acts::eBoundPhi);
      varBound[Acts::eBoundTheta] =
          fullVar(Acts::eBoundTheta, Acts::eBoundTheta);
      varBound[Acts::eBoundTime] = fullVar(Acts::eBoundTime, Acts::eBoundTime);
    }

    /// Convenience function to fill the cluster information
    ///
    /// @param c The cluster
    void fillCluster(const Cluster& c) {
      nch = static_cast<int>(c.channels.size());
      cSize[0] = static_cast<int>(c.sizeLoc0);
      cSize[1] = static_cast<int>(c.sizeLoc1);
      for (auto ch : c.channels) {
        chId[0]->push_back(static_cast<int>(ch.bin[0]));
        chId[1]->push_back(static_cast<int>(ch.bin[1]));
        chValue->push_back(static_cast<float>(ch.activation));
      }
    }
  };

  /// Constructor with
  /// @param cfg configuration struct
  /// @param output logging level
  RootMeasurementWriter(const Config& config, Acts::Logging::Level level);

  /// Virtual destructor
  ~RootMeasurementWriter() override;

  /// End-of-run hook
  ProcessCode finalize() override;

  /// Get const access to the config
  const Config& config() const { return m_cfg; }

 protected:
  /// This implementation holds the actual writing method
  /// and is called by the WriterT<>::write interface
  ///
  /// @param ctx The Algorithm context with per event information
  /// @param measurements is the data to be written out
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const MeasurementContainer& measurements) override;

 private:
  Config m_cfg;
  std::mutex m_writeMutex;  ///< protect multi-threaded writes
  TFile* m_outputFile;      ///< the output file
  Acts::GeometryHierarchyMap<std::unique_ptr<DigitizationTree>>
      m_outputTrees;  ///< the output trees
  std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*>
      m_dSurfaces;  ///< All surfaces that could carry measurements

  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};
  ReadDataHandle<IndexMultimap<Index>> m_inputMeasurementSimHitsMap{
      this, "InputMeasurementSimHitsMap"};
  ReadDataHandle<ClusterContainer> m_inputClusters{this, "InputClusters"};
};

}  // namespace ActsExamples

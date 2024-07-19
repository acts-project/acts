// This file is part of the Acts project.
//
// Copyright (C) 2017-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootMeasurementWriter.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "ActsExamples/EventData/AverageSimHits.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <cstddef>
#include <ios>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <variant>

#include <TFile.h>

namespace Acts {
class Surface;
}  // namespace Acts

namespace ActsExamples {

struct RootMeasurementWriter::DigitizationTree {
  const std::array<std::string, Acts::eBoundSize> bNames = {
      "loc0", "loc1", "phi", "theta", "qop", "time"};

  TTree* tree = nullptr;

  // Identification parameters
  int eventNr = 0;
  int volumeID = 0;
  int layerID = 0;
  int surfaceID = 0;

  // Truth parameters
  float trueBound[Acts::eBoundSize] = {};
  float trueGx = 0.;
  float trueGy = 0.;
  float trueGz = 0.;
  float incidentPhi = 0.;
  float incidentTheta = 0.;

  // Reconstruction information
  float recBound[Acts::eBoundSize] = {};
  float varBound[Acts::eBoundSize] = {};

  // Cluster information comprised of
  // nch :  number of channels
  // cSize : cluster size in loc0 and loc1
  // chId : channel identification
  // chValue: value/activation of the channel
  int nch = 0;
  int cSize[2] = {};
  std::array<std::vector<int>, 2> chId;
  std::vector<float> chValue;

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
    for (unsigned int ib = 0; ib < Acts::eBoundSize; ++ib) {
      tree->Branch(("true_" + bNames[ib]).c_str(), &trueBound[ib]);
    }
    tree->Branch("true_x", &trueGx);
    tree->Branch("true_y", &trueGy);
    tree->Branch("true_z", &trueGz);
    tree->Branch("true_incident_phi", &incidentPhi);
    tree->Branch("true_incident_theta", &incidentTheta);
  }

  /// Constructor from tree name
  DigitizationTree(const std::string& treeName) { setupTree(treeName); }

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
    tree->Branch("clus_size", &nch);
    tree->Branch("channel_value", &chValue);
    // Both are allocated, but only relevant ones are set
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
    trueBound[Acts::eBoundQOverP] = std::numeric_limits<float>::quiet_NaN();
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
    for (unsigned int ib = 0; ib < Acts::eBoundSize; ++ib) {
      recBound[ib] = fullVect[ib];
    }

    Acts::BoundSquareMatrix fullVar =
        m.expander() * m.covariance() * m.expander().transpose();
    for (unsigned int ib = 0; ib < Acts::eBoundSize; ++ib) {
      varBound[ib] = fullVar(ib, ib);
    }
  }

  /// Convenience function to fill the cluster information
  ///
  /// @param c The cluster
  void fillCluster(const Cluster& c) {
    nch = static_cast<int>(c.channels.size());
    cSize[0] = static_cast<int>(c.sizeLoc0);
    cSize[1] = static_cast<int>(c.sizeLoc1);
    for (auto ch : c.channels) {
      chId[0].push_back(static_cast<int>(ch.bin[0]));
      chId[1].push_back(static_cast<int>(ch.bin[1]));
      chValue.push_back(static_cast<float>(ch.activation));
    }
  }

  /// Clear the tree
  void clear() {
    for (unsigned int ib = 0; ib < Acts::eBoundSize; ++ib) {
      trueBound[ib] = 0.;
      recBound[ib] = 0.;
      varBound[ib] = 0.;
    }
    trueGx = 0.;
    trueGy = 0.;
    trueGz = 0.;
    incidentPhi = 0.;
    incidentTheta = 0.;
    nch = 0;
    cSize[0] = 0;
    cSize[1] = 0;
    chId[0].clear();
    chId[1].clear();
    chValue.clear();
  }
};

RootMeasurementWriter::RootMeasurementWriter(
    const RootMeasurementWriter::Config& config, Acts::Logging::Level level)
    : WriterT(config.inputMeasurements, "RootMeasurementWriter", level),
      m_cfg(config) {
  // Input container for measurements is already checked by base constructor
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  if (m_cfg.inputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-to-simulated-hits map input collection");
  }

  m_inputClusters.maybeInitialize(m_cfg.inputClusters);
  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputMeasurementSimHitsMap.initialize(m_cfg.inputMeasurementSimHitsMap);

  if (m_cfg.surfaceByIdentifier.empty()) {
    throw std::invalid_argument("Missing Surface-GeoID association map");
  }
  // Setup ROOT File
  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }

  m_outputFile->cd();

  m_outputTree = std::make_unique<DigitizationTree>("measurements");
}

RootMeasurementWriter::~RootMeasurementWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ProcessCode RootMeasurementWriter::finalize() {
  /// Close the file if it's yours
  m_outputFile->cd();
  m_outputTree->tree->Write();
  m_outputFile->Close();

  return ProcessCode::SUCCESS;
}

ProcessCode RootMeasurementWriter::writeT(
    const AlgorithmContext& ctx, const MeasurementContainer& measurements) {
  const auto& simHits = m_inputSimHits(ctx);
  const auto& hitSimHitsMap = m_inputMeasurementSimHitsMap(ctx);

  const ClusterContainer* clusters = nullptr;
  if (!m_cfg.inputClusters.empty()) {
    clusters = &m_inputClusters(ctx);
  }

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  for (Index hitIdx = 0u; hitIdx < measurements.size(); ++hitIdx) {
    const auto& meas = measurements[hitIdx];

    std::visit(
        [&](const auto& m) {
          Acts::GeometryIdentifier geoId =
              m.sourceLink().template get<IndexSourceLink>().geometryId();
          // find the corresponding surface
          auto surfaceItr = m_cfg.surfaceByIdentifier.find(geoId);
          if (surfaceItr == m_cfg.surfaceByIdentifier.end()) {
            return;
          }
          const Acts::Surface& surface = *(surfaceItr->second);

          // Fill the identification
          m_outputTree->fillIdentification(ctx.eventNumber, geoId);

          // Find the contributing simulated hits
          auto indices = makeRange(hitSimHitsMap.equal_range(hitIdx));
          // Use average truth in the case of multiple contributing sim hits
          auto [local, pos4, dir] = averageSimHits(ctx.geoContext, surface,
                                                   simHits, indices, logger());
          Acts::RotationMatrix3 rot =
              surface
                  .referenceFrame(ctx.geoContext, pos4.segment<3>(Acts::ePos0),
                                  dir)
                  .inverse();
          std::pair<double, double> angles =
              Acts::VectorHelpers::incidentAngles(dir, rot);

          m_outputTree->fillTruthParameters(local, pos4, dir, angles);
          m_outputTree->fillBoundMeasurement(m);
          if (clusters != nullptr) {
            const auto& c = (*clusters)[hitIdx];
            m_outputTree->fillCluster(c);
          }

          m_outputTree->tree->Fill();
          m_outputTree->clear();
        },
        meas);
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

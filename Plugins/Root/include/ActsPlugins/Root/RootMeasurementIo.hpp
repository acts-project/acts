// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "ActsPlugins/Root/detail/RootBranchPtr.hpp"

#include <array>
#include <tuple>
#include <vector>

class TTree;

namespace Acts {
class GeometryIdentifier;
}

namespace ActsPlugins {
/// @addtogroup root_plugin
/// @{

/// @brief Helper class to manage the I/O of measurements and associated clusters
/// to and from ROOT files.
class RootMeasurementIo {
 public:
  /// Configuration struct for measurement I/O
  struct Config {
    /// Indicate the reconstruction indices to be stored
    std::vector<Acts::BoundIndices>& recoIndices;

    /// Indicate the cluster indices to be stored
    std::vector<Acts::BoundIndices>& clusterIndices;
  };

  /// Constructor from configuration struct
  ///
  /// @param config the configuration for the accessor
  explicit RootMeasurementIo(const Config& config);

  /// @brief sets the branch connection for writing to a file
  ///
  /// @param measurementTree the TTree to write the measurement track to
  void connectForWrite(TTree& measurementTree);

  /// @brief sets the branch connection for reading from a file
  ///
  /// @param measurementTree the TTree to read the measurement track from
  void connectForRead(TTree& measurementTree);

  /// The event number of the currently loaded entry
  ///
  /// @note Only valid after connectForRead() and TTree::GetEntry()
  int eventNr() const { return m_measurementPayload.eventNr; }

  /// The geometry identifier of the currently loaded entry
  ///
  /// @note Only valid after connectForRead() and TTree::GetEntry()
  Acts::GeometryIdentifier geometryId() const;

  /// Reconstruct the bound measurement values, variances and subspace
  /// indices of the currently loaded entry. Since only the values for the
  /// dimensions that are actually part of the measurement are filled on
  /// write (the rest stay NaN), the subspace indices are recovered from the
  /// configured @c recoIndices by checking for non-NaN entries.
  ///
  /// @note Only valid after connectForRead() and TTree::GetEntry()
  ///
  /// @return a tuple of measurement values, variances and subspace indices
  std::tuple<std::vector<double>, std::vector<double>,
             std::vector<unsigned int>>
  boundMeasurement() const;

  /// The global position of the currently loaded entry (if available)
  ///
  /// @note Only valid after connectForRead() and TTree::GetEntry()
  Acts::Vector3 globalPosition() const;

  /// The cluster channels of the currently loaded entry (if available)
  ///
  /// @note Only valid after connectForRead() and TTree::GetEntry()
  std::vector<std::tuple<int, int, float>> clusterChannels() const;

  /// The cluster size in loc0/loc1 of the currently loaded entry
  ///
  /// @note Only valid after connectForRead() and TTree::GetEntry()
  std::array<int, 2> clusterSize() const { return m_clusterPayload.clusterSize; }

  /// Convenience function to register identification
  ///
  /// @param evnt The event number
  /// @param geoId The geometry identifier of the measurement
  void fillIdentification(int evnt, const Acts::GeometryIdentifier& geoId);

  /// Convenience function to register the truth parameters
  ///
  /// @param lp The true local position
  /// @param xt The true 4D global position
  /// @param dir The true particle direction
  /// @param angles The incident angles
  void fillTruthParameters(const Acts::Vector2& lp, const Acts::Vector4& xt,
                           const Acts::Vector3& dir,
                           const std::pair<double, double> angles);

  /// Convenience function to fill bound parameters
  ///  - abstracted to be used in different contexts
  ///
  /// @param measurement The measurement parameters
  /// @param variances The measurement variances (assumed diagonally)
  /// @param subspaceIndex The subspace indices of the measurement
  void fillBoundMeasurement(const std::vector<double>& measurement,
                            const std::vector<double>& variances,
                            const std::vector<unsigned int>& subspaceIndex);

  /// Fill global information of the cluster/measurement
  ///
  /// @param pos The global position of the cluster
  void fillGlobalPosition(const Acts::Vector3& pos);

  /// Convenience function to fill the cluster information
  ///  - abstracted to be used in different contexts
  ///
  /// @param channels The channel information
  void fillCluster(const std::vector<std::tuple<int, int, float>>& channels);

  /// Clear the payload
  void clear();

 private:
  // Names of the bound parameters
  static constexpr std::array<std::string, Acts::eBoundSize> bNames = {
      "loc0", "loc1", "phi", "theta", "qop", "time"};

  /// The configuration for the accessor
  Config m_cfg;

  struct MeasurementPayload {
    // Identification parameters
    int eventNr = 0;
    int volumeID = 0;
    int layerID = 0;
    int surfaceID = 0;
    int extraID = 0;

    // Reconstructed information
    std::array<float, Acts::eBoundSize> recBound = {};
    std::array<float, Acts::eBoundSize> varBound = {};

    float recGx = 0.;
    float recGy = 0.;
    float recGz = 0.;

    // Truth parameters
    std::array<float, Acts::eBoundSize> trueBound = {};
    float trueGx = 0.;
    float trueGy = 0.;
    float trueGz = 0.;
    float incidentPhi = 0.;
    float incidentTheta = 0.;

    // Residuals and pulls
    std::array<float, Acts::eBoundSize> residual = {};
    std::array<float, Acts::eBoundSize> pull = {};
  };

  struct ClusterPayload {
    // Cluster information
    int nch = 0;
    std::array<int, 2> clusterSize = {0, 0};
    // ROOT (in particular TChain) requires the address of a pointer to bind
    // std::vector<T> branches for reading, so these are owned indirectly.
    std::array<ActsExamples::RootBranchPtr<std::vector<int>>, 2> chId;
    ActsExamples::RootBranchPtr<std::vector<float>> chValue;
  };

  MeasurementPayload m_measurementPayload;
  ClusterPayload m_clusterPayload;
};
/// @}
}  // namespace ActsPlugins

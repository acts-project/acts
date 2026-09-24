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

#include <array>
#include <cstdint>
#include <vector>

class TTree;

namespace Acts {
class GeometryIdentifier;
}

namespace ActsPlugins {
/// @addtogroup root_plugin
/// @{

/// @brief Helper class to manage the I/O of measurements to and from ROOT files.
class RootMeasurementIo {
 public:
  /// Configuration struct for measurement I/O
  struct Config {
    /// Indicate the reconstruction indices to be stored
    std::vector<Acts::BoundIndices>& recoIndices;
  };

  /// Constructor from configuration struct
  ///
  /// @param config the configuration for the accessor
  explicit RootMeasurementIo(const Config& config);

  /// @brief sets the branch connection for writing to a file
  ///
  /// @param measurementTree the TTree to write the measurement track to
  void connectForWrite(TTree& measurementTree);

  /// Convenience function to register identification
  ///
  /// @param evnt The event number
  /// @param measurementId The index of the measurement within the event
  /// @param geoId The geometry identifier of the measurement
  void fillIdentification(int evnt, std::uint64_t measurementId,
                          const Acts::GeometryIdentifier& geoId);

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
    std::uint64_t measurementID = 0;
    int volumeID = 0;
    int layerID = 0;
    int surfaceID = 0;
    int extraID = 0;

    // Reconstructed information
    std::array<float, Acts::eBoundSize> recBound = {};
    std::array<float, Acts::eBoundSize> varBound = {};

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

  MeasurementPayload m_measurementPayload;
};
/// @}
}  // namespace ActsPlugins

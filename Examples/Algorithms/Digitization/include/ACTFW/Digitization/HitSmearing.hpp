// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <string>
#include <unordered_map>

#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "ACTFW/Framework/RandomNumbers.hpp"
#include "Acts/Geometry/GeometryID.hpp"

namespace Acts {
class Surface;
class TrackingGeometry;
}  // namespace Acts

namespace FW {

/// Create fittable measurements using truth smearing.
///
/// The truth information is smeared in the local measurement frame using
/// Gaussian noise to generate a fittable measurement, i.e. a source link.
class HitSmearing final : public BareAlgorithm {
 public:
  struct Config {
    /// Input collection of simulated hits.
    std::string inputSimulatedHits;
    /// Output collection for source links with smeared measurements.
    std::string outputSourceLinks;
    /// Width of the Gaussian smearing, i.e. resolution; must be positive.
    double sigmaLoc0 = -1;
    double sigmaLoc1 = -1;
    /// Tracking geometry required to access global-to-local transforms.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// Random numbers tool.
    std::shared_ptr<const RandomNumbers> randomNumbers = nullptr;
  };

  HitSmearing(const Config& cfg, Acts::Logging::Level lvl);

  ProcessCode execute(const AlgorithmContext& ctx) const final override;

 private:
  Config m_cfg;
  /// Lookup container for hit surfaces that generate smeared hits
  std::unordered_map<Acts::GeometryID, const Acts::Surface*> m_surfaces;
};

}  // namespace FW

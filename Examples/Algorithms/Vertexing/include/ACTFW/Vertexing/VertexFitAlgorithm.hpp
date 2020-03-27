// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "Acts/Utilities/Units.hpp"

using namespace Acts::UnitLiterals;

namespace FWE {

class VertexFitAlgorithm : public FW::BareAlgorithm {
 public:
  struct Config {
    /// Input track collection
    std::string trackCollection;

    /// The magnetic field
    Acts::Vector3D bField;

    bool doConstrainedFit = false;

    /// Vertex constraint covariance matrix
    Acts::ActsSymMatrixD<3> constraintCov =
        Acts::Vector3D(3_mm, 3_mm, 10_mm).asDiagonal();
    /// Vertex constraint position
    Acts::Vector3D constraintPos = Acts::Vector3D(0, 0, 0);
  };

  /// Constructor
  VertexFitAlgorithm(const Config& cfg,
                     Acts::Logging::Level level = Acts::Logging::INFO);

  /// Framework execute method
  /// @param [in] context is the Algorithm context for event consistency
  FW::ProcessCode execute(
      const FW::AlgorithmContext& context) const final override;

 private:
  /// The config class
  Config m_cfg;
};

}  // namespace FWE

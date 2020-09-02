// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"

/// Set the identifier PLUGIN
#ifdef ACTS_CORE_IDENTIFIER_PLUGIN
#include ACTS_CORE_IDENTIFIER_PLUGIN
#else
using Identifier = Acts::MinimalSourceLink;
#endif

#include "Acts/Plugins/Digitization/DigitizationCell.hpp"
#include "Acts/Plugins/Digitization/DigitizationModule.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {

template <BoundIndices... params>
using Measurement_t = Measurement<Identifier, BoundIndices, params...>;

class PlanarModuleCluster
    : public Measurement_t<BoundIndices::eBoundLoc0, BoundIndices::eBoundLoc1,
                           BoundIndices::eBoundTime> {
 public:
  /// Constructor from DigitizationCells
  ///
  /// @param [in] mSurface is the module surface
  /// @param [in] cIdentifier is the channel identifier of the local position
  /// @param [in] cov is the covariance matrix
  /// @param [in] loc0 is the local position in the first coordinate
  /// @param [in] loc1 is the local position in the second coordinate
  /// @param [in] t Timestamp of the cluster
  /// @param [in] dCells is the vector of digitization cells
  PlanarModuleCluster(std::shared_ptr<const Surface> mSurface,
                      const Identifier& identifier, ActsSymMatrixD<3> cov,
                      double loc0, double loc1, double t,
                      std::vector<DigitizationCell> dCells,
                      const DigitizationModule* dModule = nullptr)
      : Measurement_t<BoundIndices::eBoundLoc0, BoundIndices::eBoundLoc1,
                      BoundIndices::eBoundTime>(
            std::move(mSurface),
            identifier,  // original measurement
            std::move(cov), loc0, loc1, t),
        m_digitizationCells(std::move(dCells)),
        m_digitizationModule(dModule) {}

  /// access to the digitization cells
  ///
  /// @return the vector to the digitization cells
  const std::vector<DigitizationCell>& digitizationCells() const;

  /// access to the digitization module
  ///
  /// @return the pointer to the digitization module
  const DigitizationModule* digitizationModule() const;

 private:
  std::vector<DigitizationCell> m_digitizationCells;  /// the digitization cells
  const DigitizationModule* m_digitizationModule;  /// the digitization module
};

inline const std::vector<DigitizationCell>&
PlanarModuleCluster::digitizationCells() const {
  return m_digitizationCells;
}

inline const DigitizationModule* PlanarModuleCluster::digitizationModule()
    const {
  return m_digitizationModule;
}
}  // namespace Acts

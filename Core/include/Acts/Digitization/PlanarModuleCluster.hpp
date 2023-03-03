// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Digitization/DigitizationCell.hpp"
#include "Acts/Digitization/DigitizationModule.hpp"
#include "Acts/Digitization/DigitizationSourceLink.hpp"
#include "Acts/EventData/Measurement.hpp"

#include <array>
#include <cassert>

namespace Acts {

class PlanarModuleCluster : public Measurement<BoundIndices, 3> {
  using Base = Measurement<BoundIndices, 3>;

  static constexpr std::array<BoundIndices, 3> kIndices = {
      eBoundLoc0, eBoundLoc1, eBoundTime};

 public:
  /// Constructor from DigitizationCells
  ///
  /// @param [in] surface The surface the cluster is on
  /// @param [in] sourceLink is the link to the truth information
  /// @param [in] cov is the covariance matrix
  /// @param [in] loc0 is the local position in the first coordinate
  /// @param [in] loc1 is the local position in the second coordinate
  /// @param [in] t Timestamp of the cluster
  /// @param [in] dCells is the vector of digitization cells
  /// @param [in] dModule an optional pointer to a digitization configuration
  PlanarModuleCluster(std::shared_ptr<const Surface> surface,
                      SourceLink sourceLink, const Base::CovarianceMatrix& cov,
                      double loc0, double loc1, double t,
                      std::vector<DigitizationCell> dCells,
                      const DigitizationModule* dModule = nullptr)
      : Base(std::move(sourceLink), kIndices,
             Base::ParametersVector(loc0, loc1, t), cov),
        m_surface(std::move(surface)),
        m_digitizationCells(std::move(dCells)),
        m_digitizationModule(dModule) {
    assert(m_surface);
  }

  /// Module surface.
  const Surface& referenceObject() const { return *m_surface; }

  /// access to the digitization cells
  ///
  /// @return the vector to the digitization cells
  const std::vector<DigitizationCell>& digitizationCells() const;

  /// access to the digitization module
  ///
  /// @return the pointer to the digitization module
  const DigitizationModule* digitizationModule() const;

 private:
  std::shared_ptr<const Surface> m_surface;
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

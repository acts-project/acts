// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/SimHit.hpp"

namespace ActsExamples {

/// @brief DigitizedHit class that acts as a source link
/// in the `FittableMeasurement` for further usage
///
/// This class thus needs to comply with the `SourceLinkConcept`
class DigitizedHit {
 public:
  /// Constructor from surface
  ///
  /// @param surface is the Surface where the digitization result
  ///        has been expressed wrt
  DigitizedHit(const Acts::Surface& surface, std::vector<SimHit> simhits = {});
  DigitizedHit() = default;
  DigitizedHit(const DigitizedHit& other) = default;
  virtual ~DigitizedHit() = default;

  /// Equality operator is required by the SourceLinkConcept
  bool operator==(const DigitizedHit& other) const;

  /// Reference Surface required by the SourceLinkConcept
  const Acts::Surface& referenceSurface() const;

  /// The simulated hits - @todo change to index to container?
  const std::vector<SimHit>& simulatedHits() const;

 private:
  std::shared_ptr<const Acts::Surface> m_surface = nullptr;
  std::vector<SimHit> m_simulatedHits = {};
};

}  // namespace ActsExamples
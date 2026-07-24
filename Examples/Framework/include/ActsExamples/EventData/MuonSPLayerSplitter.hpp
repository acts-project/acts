// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Definitions/Common.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"

#include <iostream>
#include <memory>

namespace ActsExamples {
/** @brief Helper class to sort the MuonspacePoints whether they're Mdt hits, i.e. straw hits or whether they're ordinary strip hits.\
 *        Hits in each category are additionally sorted by straw or tube layer.
 * Each layer is represented by a dedicated hit vector. The sorter then holds
 * the sorted vectors in two jagged vectors and returns them via its interface.
 */
class MuonSPLayerSplitter {
 public:
  /** @brief abbreviation of the space point vector */
  using SpVec_t = std::vector<const MuonSpacePoint*>;
  /** @brief abbreviation of the layer vector. Each element inside the vector represents
   *         all hits in a given detector layer */
  using LayerVec = std::vector<SpVec_t>;

  /** @brief Constructor taking a list of hits in the same muon station. */
  MuonSPLayerSplitter(const SpVec_t& spacePoints);
  /** @brief Returns all straw hits sorted by straw layer */
  const LayerVec& strawHits() const { return m_strawHits; }
  /** @brief Returns all strip hits sorted by strip layer */
  const LayerVec& stripHits() const { return m_stripHits; }

 private:
  LayerVec m_strawHits{};
  LayerVec m_stripHits{};
};

}  // namespace ActsExamples

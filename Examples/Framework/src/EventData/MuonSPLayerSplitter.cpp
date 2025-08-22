// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#include "ActsExamples/EventData/MuonSPLayerSplitter.hpp"

namespace {
template <typename T>
void removeEmptyLayers(std::vector<std::vector<T>>& vec) {
  vec.erase(
      std::remove_if(vec.begin(), vec.end(),
                     [](const std::vector<T>& lay) { return lay.empty(); }),
      vec.end());
}
}  // namespace

namespace ActsExamples {

MuonSPLayerSplitter::MuonSPLayerSplitter(const SpVec_t& spacePoints) {
  m_strawHits.reserve(spacePoints.size());
  m_stripHits.reserve(spacePoints.size());
  using MuonId = MuonSpacePoint::MuonId;
  for (const MuonSpacePoint* sp : spacePoints) {
    LayerVec& pushMe{sp->id().technology() == MuonId::TechField::Mdt
                         ? m_strawHits
                         : m_stripHits};
    assert(sp->id().detLayer() >= 1);
    const unsigned idx = sp->id().detLayer() - 1;
    assert(idx >= 0);
    /// Ensure that there's enough space
    if (idx >= pushMe.size()) {
      pushMe.resize(idx + 1);
    }
    pushMe[idx].push_back(sp);
  }

  /** The gas gap number is not necessarily continuous. E.g, the bucket
   *         may have hits in layer 1, 2, 4 but not in layer 3. */
  removeEmptyLayers(m_strawHits);
  removeEmptyLayers(m_stripHits);
}

}  // namespace ActsExamples

// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/AccumulatedMaterialProperties.hpp"

#include "Acts/Material/detail/AverageMaterials.hpp"

void Acts::AccumulatedMaterialProperties::accumulate(MaterialProperties slab,
                                                     float pathCorrection) {
  // scale the recorded material to the equivalence contribution along the
  // surface normal
  slab.scaleThickness(1 / pathCorrection);
  m_trackAverage = detail::combineSlabs(m_trackAverage, slab);
}

void Acts::AccumulatedMaterialProperties::trackAverage(bool useEmptyTrack) {
  // average only real tracks or if empty tracks are allowed.
  if (useEmptyTrack or (0 < m_trackAverage.thickness())) {
    if (m_totalCount == 0u) {
      m_totalAverage = m_trackAverage;
    } else {
      double weightTotal = m_totalCount / (m_totalCount + 1.0);
      double weightTrack = 1 / (m_totalCount + 1.0);
      // average such that each track contributes equally.
      MaterialProperties fromTotal(m_totalAverage.material(),
                                   weightTotal * m_totalAverage.thickness());
      MaterialProperties fromTrack(m_trackAverage.material(),
                                   weightTrack * m_trackAverage.thickness());
      m_totalAverage = detail::combineSlabs(fromTotal, fromTrack);
    }
    m_totalCount += 1;
  }
  // reset track average
  m_trackAverage = MaterialProperties();
}

std::pair<Acts::MaterialProperties, unsigned int>
Acts::AccumulatedMaterialProperties::totalAverage() const {
  return {m_totalAverage, m_totalCount};
}

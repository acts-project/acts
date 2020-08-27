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
  // TODO the previous manual calculations used the regular thickness as the
  //   per-slab weight but then used the corrected total path as the
  //   normalization. this is inconsistent and i am not sure what the right
  //   approach is here?
  slab.scaleThickness(pathCorrection);
  m_trackAverage = detail::combineSlabs(m_trackAverage, slab);
}

void Acts::AccumulatedMaterialProperties::trackAverage(bool useEmptyTrack) {
  // average only real tracks or if empty tracks are allowed.
  if (useEmptyTrack or (0 < m_trackAverage.thickness())) {
    // average such that each track contributes equally.
    MaterialProperties weightedTotal(m_totalAverage.material(), m_totalCount);
    MaterialProperties weightedTrack(m_trackAverage.material(), 1);
    m_totalAverage = detail::combineSlabs(weightedTotal, weightedTrack);
    m_totalCount += 1;
  }
  // reset track average
  m_trackAverage = MaterialProperties();
}

std::pair<Acts::MaterialProperties, unsigned int>
Acts::AccumulatedMaterialProperties::totalAverage() const {
  // use averaged material properties but for unit thickness
  return {MaterialProperties(m_totalAverage.material(), 1), m_totalCount};
}

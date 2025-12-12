// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackStateProxy.hpp"

namespace Acts {

template <typename D, std::size_t M, bool ReadOnly>
inline TrackStateProxy<D, M, ReadOnly>::TrackStateProxy(
    const_if_t<ReadOnly, MultiTrajectory<D>>& trajectory, IndexType istate)
    : m_traj(&trajectory), m_istate(istate) {}

template <typename D, std::size_t M, bool ReadOnly>
inline auto TrackStateProxy<D, M, ReadOnly>::getUncalibratedSourceLink() const
    -> SourceLink {
  assert(has<hashString("uncalibratedSourceLink")>());
  return m_traj->getUncalibratedSourceLink(m_istate);
}

}  // namespace Acts

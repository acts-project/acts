// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/KalmanVertexUpdater.hpp"
#include "Acts/Vertexing/detail/KalmanVertexUpdaterImpl.hpp"

namespace Acts::KalmanVertexUpdater::detail {

template void updateVertexWithTrackImpl<4>(Vertex& vtx, TrackAtVertex& trk,
                                           int sign);

template void updateTrackWithVertexImpl<4>(TrackAtVertex& track,
                                           const Vertex& vtx);

}  // namespace Acts::KalmanVertexUpdater::detail

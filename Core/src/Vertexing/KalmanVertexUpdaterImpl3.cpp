// SPDX-PackageName: "ACTS"
// SPDX-FileCopyrightText: 2016 CERN
// SPDX-License-Identifier: MPL-2.0

#include "Acts/Vertexing/KalmanVertexUpdater.hpp"
#include "Acts/Vertexing/detail/KalmanVertexUpdaterImpl.hpp"

template void Acts::KalmanVertexUpdater::detail::updateVertexWithTrackImpl<3>(
    Vertex& vtx, TrackAtVertex& trk, int sign);

template void Acts::KalmanVertexUpdater::detail::updateTrackWithVertexImpl<3>(
    TrackAtVertex& track, const Vertex& vtx);

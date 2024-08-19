// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/KalmanVertexUpdater.hpp"
#include "Acts/Vertexing/detail/KalmanVertexUpdaterImpl.hpp"

template void Acts::KalmanVertexUpdater::detail::updateVertexWithTrackImpl<3>(
    Vertex& vtx, TrackAtVertex& trk, int sign);

template void Acts::KalmanVertexUpdater::detail::updateTrackWithVertexImpl<3>(
    TrackAtVertex& track, const Vertex& vtx);

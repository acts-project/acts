// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/TrackFindingExaTrkX/TrackFindingAlgorithmExaTrkX.hpp"

namespace ActsExamples {

std::vector<float> createFeatures(
    const SimSpacePointContainer &spacepoints, const ClusterContainer *clusters,
    const std::vector<TrackFindingAlgorithmExaTrkX::NodeFeature> &nodeFeatures,
    const std::vector<float> &featureScales);

}

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/SpacePoint.hpp"
#include "ActsExamples/TrackFindingGnn/TrackFindingAlgorithmGnn.hpp"

namespace ActsExamples {

std::vector<float> createFeatures(
    const ConstSpacePointSubset &spacePoints, const ClusterContainer *clusters,
    const std::vector<TrackFindingAlgorithmGnn::NodeFeature> &nodeFeatures,
    const std::vector<float> &featureScales);

}

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/EventData/Jets.hpp"

#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"

namespace ActsExamples {
// Initiate a truth jet with 4-momentum
const Acts::Vector4 testVec4(10.0, 0.0, 10.0, 14.1421);
const TruthJet testJet(testVec4, JetLabel::Unknown);
}  // namespace ActsExamples

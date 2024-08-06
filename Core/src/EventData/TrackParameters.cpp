// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/TrackParameters.hpp"

#include "Acts/EventData/TrackParametersConcept.hpp"

namespace Acts {

// ensure concrete classes satisfy the concepts

static_assert(
    Concepts::BoundTrackParameters<SinglyChargedBoundTrackParameters>);
static_assert(
    Concepts::BoundTrackParameters<SinglyChargedCurvilinearTrackParameters>);
static_assert(Concepts::FreeTrackParameters<SinglyChargedFreeTrackParameters>);

static_assert(Concepts::BoundTrackParameters<NeutralBoundTrackParameters>);
static_assert(
    Concepts::BoundTrackParameters<NeutralCurvilinearTrackParameters>);
static_assert(Concepts::FreeTrackParameters<NeutralFreeTrackParameters>);

static_assert(Concepts::BoundTrackParameters<BoundTrackParameters>);
static_assert(Concepts::BoundTrackParameters<CurvilinearTrackParameters>);
static_assert(Concepts::FreeTrackParameters<FreeTrackParameters>);

}  // namespace Acts

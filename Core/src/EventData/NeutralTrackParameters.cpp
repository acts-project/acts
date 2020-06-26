// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/NeutralTrackParameters.hpp"

#include "Acts/EventData/TrackParametersConcept.hpp"

namespace Acts {

// explicitly instantiate templates
template class SingleBoundTrackParameters<Neutral>;
template class SingleCurvilinearTrackParameters<Neutral>;
template class SingleFreeTrackParameters<Neutral>;

// ensure concrete classes satisfy the concepts
static_assert(
    Concepts::BoundTrackParametersConcept<NeutralBoundTrackParameters>);
static_assert(
    Concepts::BoundTrackParametersConcept<NeutralCurvilinearTrackParameters>);
static_assert(Concepts::FreeTrackParametersConcept<NeutralFreeTrackParameters>);

}  // namespace Acts

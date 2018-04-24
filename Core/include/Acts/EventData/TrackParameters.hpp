// This file is part of the ACTS project.
//
// Copyright (C) 2016-2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_TRACKPARAMETERS_H
#define ACTS_TRACKPARAMETERS_H

#include "Acts/EventData/ChargePolicy.hpp"
#include "Acts/EventData/SingleBoundTrackParameters.hpp"
#include "Acts/EventData/SingleCurvilinearTrackParameters.hpp"
#include "Acts/EventData/SingleTrackParameters.hpp"

namespace Acts {
typedef SingleTrackParameters<ChargedPolicy>            TrackParameters;
typedef SingleCurvilinearTrackParameters<ChargedPolicy> CurvilinearParameters;
typedef SingleBoundTrackParameters<ChargedPolicy>       BoundParameters;
}  // end of namespace Acts

#endif  // ACTS_TRACKPARAMETERS_H

// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_TRACKPARAMETERS_H
#define ACTS_TRACKPARAMETERS_H

#include "ACTS/EventData/ChargePolicy.hpp"
#include "ACTS/EventData/SingleBoundTrackParameters.hpp"
#include "ACTS/EventData/SingleCurvilinearTrackParameters.hpp"
#include "ACTS/EventData/SingleTrackParameters.hpp"

namespace Acts {
typedef SingleTrackParameters<ChargedPolicy>            TrackParameters;
typedef SingleCurvilinearTrackParameters<ChargedPolicy> CurvilinearParameters;
typedef SingleBoundTrackParameters<ChargedPolicy>       BoundParameters;
}  // end of namespace Acts

#endif  // ACTS_TRACKPARAMETERS_H

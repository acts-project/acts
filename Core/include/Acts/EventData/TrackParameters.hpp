// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Charge.hpp"
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/GenericCurvilinearTrackParameters.hpp"
#include "Acts/EventData/GenericFreeTrackParameters.hpp"

namespace Acts {

using BoundTrackParameters = GenericBoundTrackParameters<SinglyCharged>;
using CurvilinearTrackParameters =
    GenericCurvilinearTrackParameters<SinglyCharged>;
using FreeTrackParameters = GenericFreeTrackParameters<SinglyCharged>;

using NeutralBoundTrackParameters = GenericBoundTrackParameters<Neutral>;
using NeutralCurvilinearTrackParameters =
    GenericCurvilinearTrackParameters<Neutral>;
using NeutralFreeTrackParameters = GenericFreeTrackParameters<Neutral>;

}  // namespace Acts

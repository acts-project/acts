// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Charge.hpp"
#include "Acts/EventData/SingleBoundTrackParameters.hpp"
#include "Acts/EventData/SingleCurvilinearTrackParameters.hpp"
#include "Acts/EventData/SingleFreeTrackParameters.hpp"

namespace Acts {

extern template class SingleBoundTrackParameters<Neutral>;
extern template class SingleCurvilinearTrackParameters<Neutral>;
extern template class SingleFreeTrackParameters<Neutral>;

using NeutralBoundTrackParameters = SingleBoundTrackParameters<Neutral>;
using NeutralCurvilinearTrackParameters =
    SingleCurvilinearTrackParameters<Neutral>;
using NeutralFreeTrackParameters = SingleFreeTrackParameters<Neutral>;

}  // namespace Acts

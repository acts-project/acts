// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/ChargePolicy.hpp"
#include "Acts/EventData/SingleBoundTrackParameters.hpp"
#include "Acts/EventData/SingleCurvilinearTrackParameters.hpp"
#include "Acts/EventData/SingleFreeTrackParameters.hpp"

namespace Acts {

extern template class SingleBoundTrackParameters<NeutralPolicy>;
extern template class SingleCurvilinearTrackParameters<NeutralPolicy>;
extern template class SingleFreeTrackParameters<NeutralPolicy>;

using NeutralBoundTrackParameters = SingleBoundTrackParameters<NeutralPolicy>;
using NeutralCurvilinearTrackParameters =
    SingleCurvilinearTrackParameters<NeutralPolicy>;
using NeutralFreeTrackParameters = SingleFreeTrackParameters<NeutralPolicy>;

}  // namespace Acts

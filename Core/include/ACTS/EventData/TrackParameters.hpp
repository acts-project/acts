// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_TRACKPARAMETERS_H
#define ACTS_TRACKPARAMETERS_H 1

#include "ACTS/EventData/SingleCurvilinearTrackParameters.hpp"
#include "ACTS/EventData/SingleBoundTrackParameters.hpp"
#include "ACTS/EventData/SingleTrackParameters.hpp"
#include "ACTS/EventData/ChargePolicy.hpp"

namespace Acts
{
    typedef SingleTrackParameters<ChargedPolicy> TrackParameters;
    typedef SingleCurvilinearTrackParameters<ChargedPolicy> CurvilinearParameters;
    typedef SingleBoundTrackParameters<ChargedPolicy> BoundParameters;
} // end of namespace Acts


/**Overload of << operator for both, MsgStream and std::ostream for debug output*/
std::ostream& operator << ( std::ostream& sl, const Acts::TrackParameters& pars);

#endif // ACTS_TRACKPARAMETERS_H

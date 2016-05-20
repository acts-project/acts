// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/EventData/TrackParameters.hpp"

namespace Acts
{
 template class SingleTrackParameters<ChargedPolicy>;
 template class SingleCurvilinearTrackParameters<ChargedPolicy>;
 template class SingleBoundTrackParameters<ChargedPolicy>;
} // end of namespace Acts

/**Overload of << operator for both, MsgStream and std::ostream for debug output*/
std::ostream& operator << ( std::ostream& sl, const Acts::TrackParameters& pars)
{ return pars.dump(sl); }

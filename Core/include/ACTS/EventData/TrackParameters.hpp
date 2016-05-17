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

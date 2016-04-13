#ifndef ACTS_TRACKPARAMETERS_H
#define ACTS_TRACKPARAMETERS_H 1

#include "EventData/SingleCurvilinearTrackParameters.h"
#include "EventData/SingleBoundTrackParameters.h"
#include "EventData/SingleTrackParameters.h"
#include "EventData/ChargePolicy.h"

namespace Acts
{
    typedef SingleTrackParameters<ChargedPolicy> TrackParameters;
    typedef SingleCurvilinearTrackParameters<ChargedPolicy> CurvilinearParameters;
    typedef SingleBoundTrackParameters<ChargedPolicy> BoundParameters;
} // end of namespace Acts


/**Overload of << operator for both, MsgStream and std::ostream for debug output*/
std::ostream& operator << ( std::ostream& sl, const Acts::TrackParameters& pars);

#endif // ACTS_TRACKPARAMETERS_H

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

#include "ACTS/NeutralParameters.h"

namespace Acts
{
 template class SingleTrackParameters<NeutralPolicy>;
 template class SingleCurvilinearTrackParameters<NeutralPolicy>;
 template class SingleBoundTrackParameters<NeutralPolicy>;
} // end of namespace Acts

/**Overload of << operator for both, MsgStream and std::ostream for debug output*/
std::ostream& operator << ( std::ostream& sl, const Acts::NeutralParameters& pars)
{ return pars.dump(sl); }


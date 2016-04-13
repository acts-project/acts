#ifndef ACTS_NEUTRALPARAMETERS_H
#define ACTS_NEUTRALPARAMETERS_H 1

#include "ParametersBase/SingleTrackParameters.h"
#include "ParametersBase/SingleCurvilinearTrackParameters.h"
#include "ParametersBase/SingleBoundTrackParameters.h"
#include "ParametersBase/ChargePolicy.h"

namespace Acts
{
  typedef SingleTrackParameters<NeutralPolicy> NeutralParameters;
  typedef SingleCurvilinearTrackParameters<NeutralPolicy> NeutralCurvilinearParameters;
  typedef SingleBoundTrackParameters<NeutralPolicy> NeutralBoundParameters;
} // end of namespace Acts

/**Overload of << operator for both, MsgStream and std::ostream for debug output*/
std::ostream& operator << ( std::ostream& sl, const Acts::NeutralParameters& pars);

#endif // ACTS_NEUTRALPARAMETERS_H

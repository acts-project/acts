#include "Acts/Plugins/ExaTrkX/ExaTrkXTrackFindingBase.hpp"

namespace Acts {

void ExaTrkXTrackFindingBase::getTracks(
    std::vector<float>& inputValues, std::vector<int>& spacepointIDs,
    std::vector<std::vector<int> >& trackCandidates) const {
  auto timeInfo = ExaTrkXTime{};
  getTracks(inputValues, spacepointIDs, trackCandidates, timeInfo);
}

}  // namespace Acts

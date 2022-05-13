#pragma once

#include "ActsFatras/EventData/Hit.hpp"

#include "edm4hep/SimTrackerHit.h"

namespace Acts {

ActsFatras::Hit convertEDM4hepSimHit(const edm4hep::SimTrackerHit& sth);

} // namespace Acts

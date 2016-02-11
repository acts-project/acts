#include "GaudiKernel/DeclareFactoryEntries.h"
#include "RungeKuttaEngine/RungeKuttaEngine.h"

using namespace Ats;

DECLARE_SERVICE_FACTORY( RungeKuttaEngine )

/** factory entries need to have the name of the package */
DECLARE_FACTORY_ENTRIES( RungeKuttaEngine )
{
    DECLARE_SERVICE( RungeKuttaEngine )
}


#include "GaudiKernel/DeclareFactoryEntries.h"
#include "TrkExRungeKuttaEngine/RungeKuttaEngine.h"

using namespace Ats;

DECLARE_SERVICE_FACTORY( RungeKuttaEngine )

/** factory entries need to have the name of the package */
DECLARE_FACTORY_ENTRIES( TrkExRungeKuttaEngine )
{
    DECLARE_SERVICE( RungeKuttaEngine )
}


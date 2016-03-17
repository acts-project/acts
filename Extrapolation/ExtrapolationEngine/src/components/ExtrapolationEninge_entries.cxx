#include "GaudiKernel/DeclareFactoryEntries.h"
#include "ExtrapolationEngine/ExtrapolationEngine.h"
#include "ExtrapolationEngine/StaticEngine.h"
#include "ExtrapolationEngine/StaticNavigationEngine.h"
#include "ExtrapolationEngine/MaterialEffectsEngine.h"
#include "ExtrapolationEngine/PropagationEngine.h"

using namespace Acts;

DECLARE_SERVICE_FACTORY( ExtrapolationEngine )
DECLARE_SERVICE_FACTORY( MaterialEffectsEngine )
DECLARE_SERVICE_FACTORY( StaticEngine )
DECLARE_SERVICE_FACTORY( StaticNavigationEngine )
DECLARE_SERVICE_FACTORY( PropagationEngine )

/** factory entries need to have the name of the package */
DECLARE_FACTORY_ENTRIES( ExtrapolationEngine )
{
    DECLARE_SERVICE( ExtrapolationEngine )
    DECLARE_SERVICE( MaterialEffectsEngine )        
    DECLARE_SERVICE( StaticEngine )
    DECLARE_SERVICE( StaticNavigationEngine )
    DECLARE_SERVICE( PropagationEngine )
}

#include "GaudiKernel/DeclareFactoryEntries.h"
#include "FatrasServices/FatrasMaterialEffectsEngine.h"

using namespace Acts;

DECLARE_SERVICE_FACTORY( FatrasMaterialEffectsEngine )

/** factory entries need to have the name of the package */

DECLARE_FACTORY_ENTRIES( FatrasServices )
{
    DECLARE_SERVICE( FatrasMaterialEffectsEngine )
}

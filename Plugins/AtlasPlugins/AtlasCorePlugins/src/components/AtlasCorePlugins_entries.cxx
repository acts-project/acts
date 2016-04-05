#include "GaudiKernel/DeclareFactoryEntries.h"
#include "AtlasCorePlugins/AtlasRandomNumberSvc.h"

using namespace Acts;

DECLARE_SERVICE_FACTORY( AtlasRandomNumberSvc )

/** factory entries need to have the name of the package */
DECLARE_FACTORY_ENTRIES( AtlasCorePlugins )
{
    DECLARE_SERVICE( AtlasRandomNumberSvc )
}

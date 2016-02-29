#include "GaudiKernel/DeclareFactoryEntries.h"

#include "GenericGeometryTools/GenericLayerBuilder.h"

using namespace Ats;

DECLARE_TOOL_FACTORY( GenericLayerBuilder )


/** factory entries need to have the name of the package */
DECLARE_FACTORY_ENTRIES( GenericGeometryTools )
{
  DECLARE_TOOL( GenericLayerBuilder )
    
}

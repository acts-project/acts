#include "GaudiKernel/DeclareFactoryEntries.h"
#include "JsonWriters/GeometryJsonWriter.h"
#include "JsonWriters/ParametersJsonWriter.h"


using namespace Ats;

DECLARE_TOOL_FACTORY( GeometryJsonWriter )
DECLARE_TOOL_FACTORY( ParametersJsonWriter )

/** factory entries need to have the name of the package */
DECLARE_FACTORY_ENTRIES( JsonWriters )
{
    DECLARE_TOOL( GeometryJsonWriter )
    DECLARE_TOOL( ParametersJsonWriter )
    
}


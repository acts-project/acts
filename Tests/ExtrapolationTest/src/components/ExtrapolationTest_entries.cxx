#include "GaudiKernel/DeclareFactoryEntries.h"
#include "ExtrapolationTest/ExtrapolationEngineTest.h"
#include "ExtrapolationTest/PropagationEngineTest.h"

using namespace Acts;

DECLARE_ALGORITHM_FACTORY( ExtrapolationEngineTest )
DECLARE_ALGORITHM_FACTORY( PropagationEngineTest )

/** factory entries need to have the name of the package */
DECLARE_FACTORY_ENTRIES( ExtrapolationTest )
{
    DECLARE_ALGORITHM( ExtrapolationEngineTest )
    DECLARE_ALGORITHM( PropagationEngineTest )        
}

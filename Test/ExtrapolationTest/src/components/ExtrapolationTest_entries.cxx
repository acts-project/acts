#include "GaudiKernel/DeclareFactoryEntries.h"
#include "ExtrapolationTest/ExtrapolationEngineTest.h"

using namespace Ats;

DECLARE_ALGORITHM_FACTORY( ExtrapolationEngineTest )

/** factory entries need to have the name of the package */
DECLARE_FACTORY_ENTRIES( ExtrapolationTest )
{
    DECLARE_ALGORITHM( ExtrapolationEngineTest )
}

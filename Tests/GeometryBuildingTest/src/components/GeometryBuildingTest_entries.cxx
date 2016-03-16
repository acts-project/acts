#include "GaudiKernel/DeclareFactoryEntries.h"
// Examples module
#include "GeometryBuildingTest/TrackingGeometryTest.h"

using namespace Acts;

DECLARE_ALGORITHM_FACTORY( TrackingGeometryTest )

/** factory entries need to have the name of the package */
DECLARE_FACTORY_ENTRIES( GeometryBuildingTest )
{
  DECLARE_ALGORITHM( TrackingGeometryTest )
}

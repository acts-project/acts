#ifndef ATLASPLUGINS_ATHENAINTERFACES_H
#define ATLASPLUGINS_ATHENAINTERFACES_H 1

// Athena version
#include "AthenaBaseComps/AthAlgorithm.h"
#include "AthenaBaseComps/AthAlgTool.h"
#include "AthenaBaseComps/AthService.h"
#include "AtlasCorePlugins/AthenaMessaging.h"

namespace Acts {
    typedef AthAlgorithm AlgorithmBase;
    typedef AthAlgTool   AlgToolBase;
    typedef AthService   ServiceBase;
}

#endif
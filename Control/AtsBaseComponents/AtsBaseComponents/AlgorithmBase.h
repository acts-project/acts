///////////////////////////////////////////////////////////////////
// AlgorithmBase.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef BASECOMPSINTERFACES_ALGORITHMBASE_H
#define BASECOMPSINTERFACES_ALGORITHMBASE_H 1

#include "AtsBaseComponents/BaseMacros.h"

// ATLAS build
#ifndef ATS_GAUDI_BUILD

// Athena version
#include "AthenaBaseComps/AthAlgorithm.h"
#include "AtsBaseComponents/MsgMacros.h"
#include "AtsBaseComponents/BaseMacros.h"


namespace Ats {
    typedef AthAlgorithm AlgorithmBase;
}

// GAUDI build
#else

#include "GaudiKernel/Algorithm.h"
#include "AtsBaseComponents/MsgBase.h"
#include "AtsBaseComponents/BaseMacros.h"

class ISvcLocator;

namespace Ats {
    
    /**  @class AlgorithmBase
         simply extend the Service class with a MsgBase */
    class AlgorithmBase : public ::Algorithm, public MsgBase {
        public:
            /** Constructor */
            AlgorithmBase(const std::string& name, ISvcLocator* pSvcLocator, const std::string& version=PACKAGE_VERSION) :
                ::Algorithm(name, pSvcLocator, version),
                MsgBase(msgSvc(), name) {}
               
            /** Virtual Destructor */   
      	    virtual ~AlgorithmBase() {}
    }
}

#endif

#endif // BASECOMPSINTERFACES_ALGORITHMBASE_H

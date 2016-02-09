///////////////////////////////////////////////////////////////////
// AlgorithmBase.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_CORE_ALGORITHMBASE_H
#define ATS_CORE_ALGORITHMBASE_H 1

#include "CoreInterfaces/BaseMacros.h"

// Take the interface from another sourse
#ifdef ATS_CORE_ALGORITHM_PLUGIN
#include ATS_CORE_ALGORITHM_PLUGIN
#else 

#include "GaudiKernel/Algorithm.h"
#include "CoreInterfaces/MsgBase.h"

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
    };
}

#endif //ATS_CORE_ALGORITHM_PLUGIN

#endif // ATS_CORE_ALGORITHMBASE_H

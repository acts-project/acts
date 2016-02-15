///////////////////////////////////////////////////////////////////
// AlgorithmBase.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_CORE_ALGORITHMBASE_H
#define ATS_CORE_ALGORITHMBASE_H 1

// (a) Take the interface from another sourse
#ifdef ATS_CORE_ALGORITHM_PLUGIN
#include ATS_CORE_ALGORITHM_PLUGIN
#else 

// (b) Define here
#include "CoreInterfaces/MsgBase.h"
#include "CoreInterfaces/MsgMacros.h"
#include "GaudiKernel/Algorithm.h"

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

// (c) used in both cases
#include "CoreInterfaces/BaseMacros.h"

#endif // ATS_CORE_ALGORITHMBASE_H

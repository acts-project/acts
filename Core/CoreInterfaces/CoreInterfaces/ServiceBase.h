///////////////////////////////////////////////////////////////////
// ServiceBase.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_CORE_SERVICEBASE_H
#define ACTS_CORE_SERVICEBASE_H 1

// (a) Take the interface from another sourse
#ifdef ACTS_CORE_SERVICE_PLUGIN
#include ACTS_CORE_SERVICE_PLUGIN
#else 

// (b) define here 
#include "GaudiKernel/Service.h"
#include "CoreInterfaces/MsgBase.h"
#include "CoreInterfaces/MsgMacros.h"

class ISvcLocator;

namespace Acts {
    
    /**  @class ServiceBase
         simply extend the Service class with a MsgBase */
    class ServiceBase : public ::Service, public MsgBase {
        public:
            /** Constructor */
            ServiceBase(const std::string& name, ISvcLocator* pSvcLocator ) :
                ::Service(name, pSvcLocator),
                MsgBase(msgSvc(), name)
            {}
               
            /** Virtual Destructor */   
      	    virtual ~ServiceBase() {}
    };
}

#endif //ACTS_CORE_SERVICE_PLUGIN

// (c) used in bose cases
#include "CoreInterfaces/BaseMacros.h"

#endif // ACTS_CORE_SERVICEBASE_H

///////////////////////////////////////////////////////////////////
// ServiceBase.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_CORE_SERVICEBASE_H
#define ATS_CORE_SERVICEBASE_H 1

// (a) Take the interface from another sourse
#ifdef ATS_CORE_SERVICE_PLUGIN
#include ATS_CORE_SERVICE_PLUGIN
#else 

// (b) define here 
#include "GaudiKernel/Service.h"
#include "CoreInterfaces/MsgBase.h"
#include "CoreInterfaces/MsgMacros.h"

class ISvcLocator;

namespace Ats {
    
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

#endif //ATS_CORE_SERVICE_PLUGIN

// (c) used in bose cases
#include "CoreInterfaces/BaseMacros.h"

#endif // ATS_CORE_SERVICEBASE_H

///////////////////////////////////////////////////////////////////
// ServiceBase.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATSBASECOMPONENTS_SERVICEBASE_H
#define ATSBASECOMPONENTS_SERVICEBASE_H 1

// ATHENA build
#ifndef ATS_GAUDI_BUILD

// Athena version
#include "AthenaBaseComps/AthService.h"
#include "AtsBaseComponents/MsgMacros.h"

namespace Ats {

    /** @typedef ServiceBase from AthService */
    typedef AthService ServiceBase;
}

// GAUDI build
#else

#include "GaudiKernel/Service.h"
#include "AtsBaseComponents/MsgBase.h"

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


#endif //ATS_GAUDI_BUILD


#endif // BASECOMPSINTERFACES_SERVICEBASE_H

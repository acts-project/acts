///////////////////////////////////////////////////////////////////
// AlgToolBase.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_BASECOMPSINTERFACES_ALGTOOLBASE_H
#define ATS_BASECOMPSINTERFACES_ALGTOOLBASE_H 1

// ATLAS build
#ifndef ATS_GAUDI_BUILD

// Athena version
#include "AthenaBaseComps/AthAlgTool.h"
#include "AtsBaseComponents/MsgMacros.h"
#include "AtsBaseComponents/BaseMacros.h"

namespace Ats {
    typedef AthAlgTool AlgToolBase;
}

// GAUDI build
#else 

#include "GaudiKernel/AlgTool.h"
#include "AtsBaseComponents/MsgBase.h"
#include "AtsBaseComponents/BaseMacros.h"
    
namespace Ats {
    
    /**  @class AlgToolBase
         simply extend the AlgTool class with a MsgBase */
    class AlgToolBase : public ::AlgTool, public MsgBase {
        public:
            /** Constructor */
            AlgToolBase( const std::string& type, const std::string& name, const IInterface* parent) :
                ::AlgTool(type, name, parent),
                MsgBase(msgSvc(), name)
            {}
               
            /** Virtual Destructor */   
      	    virtual ~AlgToolBase() {}
    }
}
    
#endif // ATLAS-GAUDI-build    

#endif // BASECOMPSINTERFACES_ALGTOOLBASE_H

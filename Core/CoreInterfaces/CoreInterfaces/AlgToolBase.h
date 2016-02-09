///////////////////////////////////////////////////////////////////
// AlgToolBase.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_CORE_ALGTOOLBASE_H
#define ATS_CORE_ALGTOOLBASE_H 1

#include "CoreInterfaces/BaseMacros.h"

// Take the interface from another sourse
#ifdef ATS_CORE_ALGTOOL_PLUGIN
#include ATS_CORE_ALGTOOL_PLUGIN
#else 

#include "CoreInterfaces/MsgBase.h"
#include "GaudiKernel/AlgTool.h"
    
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
    };
}
    
#endif // ATS_CORE_ALGTOOL_PLUGIN

#endif // ATS_CORE_ALGTOOLBASE_H

///////////////////////////////////////////////////////////////////
// AlgToolBase.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_CORE_ALGTOOLBASE_H
#define ATS_CORE_ALGTOOLBASE_H 1

// (a) Take the interface from another sourse
#ifdef ATS_CORE_ALGTOOL_PLUGIN
#include ATS_CORE_ALGTOOL_PLUGIN
#else 

// (b) Define here
#include "CoreInterfaces/MsgBase.h"
#include "CoreInterfaces/MsgMacros.h"
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

// (c) use in both cases
#include "CoreInterfaces/BaseMacros.h"

#endif // ATS_CORE_ALGTOOLBASE_H

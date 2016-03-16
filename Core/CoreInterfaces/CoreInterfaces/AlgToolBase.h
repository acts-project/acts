///////////////////////////////////////////////////////////////////
// AlgToolBase.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_CORE_ALGTOOLBASE_H
#define ACTS_CORE_ALGTOOLBASE_H 1

// (a) Take the interface from another sourse
#ifdef ACTS_CORE_ALGTOOL_PLUGIN
#include ACTS_CORE_ALGTOOL_PLUGIN
#else 

// (b) Define here
#include "CoreInterfaces/MsgBase.h"
#include "CoreInterfaces/MsgMacros.h"
#include "GaudiKernel/AlgTool.h"
    
namespace Acts {
    
    /**  @class AlgToolBase
         simply extend the AlgTool class with a MsgBase */
    class AlgToolBase : public ::AlgTool, public MsgBase {
        public:
            /** Constructor */
            AlgToolBase( const std::string& type, const std::string& name, const IInterface* parent) :
                ::AlgTool(type, name, parent),
                MsgBase(msgSvc(), name)
            {}
        
            virtual StatusCode initialize() override
            {
                return ::AlgTool::initialize();
            }
            /** Virtual Destructor */
      	    virtual ~AlgToolBase() {}
    };
}

#endif // ACTS_CORE_ALGTOOL_PLUGIN

// (c) use in both cases
#include "CoreInterfaces/BaseMacros.h"

#endif // ACTS_CORE_ALGTOOLBASE_H

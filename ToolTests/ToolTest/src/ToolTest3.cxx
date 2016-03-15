///////////////////////////////////////////////////////////////////
// ToolTest3.cxx, ATS project
///////////////////////////////////////////////////////////////////

// Core module
// Geometry module
#include "ToolTest/ToolTest3.h"
#include "GaudiKernel/MsgStream.h"

DECLARE_TOOL_FACTORY(ToolTest3)

// constructor
ToolTest3::ToolTest3(const std::string& t, const std::string& n, const IInterface* p) :
Ats::AlgToolBase(t,n,p)
{
    declareInterface<IToolTest3>(this);
    MSG_INFO("ToolTest3 constructor");
}

// destructor
ToolTest3::~ToolTest3()
{}

// initialize
StatusCode ToolTest3::initialize()
{
    MSG_INFO( "initialize()" );
    StatusCode sc = Ats::AlgToolBase::initialize();
    if(!sc.isSuccess()){
        MSG_FATAL("Could not initialize Tool");
        return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
}

//finalize
StatusCode ToolTest3::finalize()
{
    MSG_INFO( "finalize()" );
    return StatusCode::SUCCESS;
}


StatusCode ToolTest3::toolTest3() const
{
    MSG_INFO("toolTest3()");
    return StatusCode::SUCCESS;
}
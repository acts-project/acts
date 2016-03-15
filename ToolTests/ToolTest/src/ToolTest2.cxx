///////////////////////////////////////////////////////////////////
// ToolTest2.cxx, ATS project
///////////////////////////////////////////////////////////////////

// Core module
// Geometry module
#include "ToolTest/ToolTest2.h"
#include "GaudiKernel/MsgStream.h"

DECLARE_TOOL_FACTORY(ToolTest2)

// constructor
ToolTest2::ToolTest2(const std::string& t, const std::string& n, const IInterface* p) :
Ats::AlgToolBase(t,n,p),
m_toolTest3(this)
{
    declareInterface<IToolTest2>(this);
    MSG_INFO("ToolTest2 constructor");
    MSG_INFO("ToolTest constructor");
    declareProperty("ToolTest3",m_toolTest3);
}

// destructor
ToolTest2::~ToolTest2()
{}

// initialize
StatusCode ToolTest2::initialize()
{
    MSG_INFO( "initialize()" );
    StatusCode sc = Ats::AlgToolBase::initialize();
    if(!sc.isSuccess()){
        MSG_FATAL("Could not initialize Tool");
        return StatusCode::FAILURE;
    }
    RETRIEVE_FATAL(m_toolTest3);
    return StatusCode::SUCCESS;
}

//finalize
StatusCode ToolTest2::finalize()
{
    MSG_INFO( "finalize()" );
    return StatusCode::SUCCESS;
}


StatusCode ToolTest2::toolTest2() const
{
    MSG_INFO("toolTest2()");
    return m_toolTest3->toolTest3();
}
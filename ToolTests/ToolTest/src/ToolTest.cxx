///////////////////////////////////////////////////////////////////
// ToolTest.cxx, ATS project
///////////////////////////////////////////////////////////////////

// Core module
// Geometry module
#include "ToolTest/ToolTest.h"
#include "GaudiKernel/MsgStream.h"

DECLARE_TOOL_FACTORY(ToolTest)

// constructor
ToolTest::ToolTest(const std::string& t, const std::string& n, const IInterface* p) :
Ats::AlgToolBase(t,n,p),
m_toolTest2(this)
{
    declareInterface<IToolTest>(this);
    MSG_INFO("ToolTest constructor");
    declareProperty("ToolTest2",m_toolTest2);
}

// destructor
ToolTest::~ToolTest()
{}

// initialize
StatusCode ToolTest::initialize()
{
    MSG_INFO( "initialize()" );
    StatusCode sc = Ats::AlgToolBase::initialize();
    if(!sc.isSuccess()){
        MSG_FATAL("Could not initialize Tool");
        return StatusCode::FAILURE;
    }
    RETRIEVE_FATAL(m_toolTest2);
    return StatusCode::SUCCESS;
}

//finalize
StatusCode ToolTest::finalize()
{
    MSG_INFO( "finalize()" );
    return StatusCode::SUCCESS;
}


StatusCode ToolTest::toolTest() const
{
    MSG_INFO("toolTest()");
    return m_toolTest2->toolTest2();
}
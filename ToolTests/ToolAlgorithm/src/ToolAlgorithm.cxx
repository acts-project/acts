//////////////////////////////////////////////////////////////////
// ToolAlgorithm.cxx, ATS project
///////////////////////////////////////////////////////////////////

// Examples module
#include "ToolAlgorithm/ToolAlgorithm.h"

DECLARE_COMPONENT(ToolAlgorithm)

ToolAlgorithm::ToolAlgorithm(const std::string& name, ISvcLocator* pSvcLocator) :
 Ats::AlgorithmBase(name, pSvcLocator),
 m_toolTest(this)
 {
     declareProperty("ToolTest",m_toolTest);
 }

StatusCode ToolAlgorithm::initialize()
{
    MSG_INFO("initialize()");
    RETRIEVE_FATAL(m_toolTest);
    return StatusCode::SUCCESS;
}

StatusCode ToolAlgorithm::finalize()
{
    MSG_INFO("finalize()");
    return StatusCode::SUCCESS;
}

StatusCode ToolAlgorithm::execute()
{
    MSG_INFO("execute()");
    return m_toolTest->toolTest();
    //return StatusCode::SUCCESS;
}

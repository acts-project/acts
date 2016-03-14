///////////////////////////////////////////////////////////////////
// ToolTest.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_GEOMETRYTOOLS_TOOLTEST_H
#define ATS_GEOMETRYTOOLS_TOOLTEST_H 1

// Core module
#include "CoreInterfaces/AlgToolBase.h"
// Geoemtry module
#include "ToolTestInterfaces/IToolTest.h"
// Gaudi & Athena
#include "GaudiKernel/ToolHandle.h"
#include "ToolTestInterfaces/IToolTest2.h"

    
class ToolTest : public Ats::AlgToolBase, virtual public IToolTest {
        
    public:
        /** constructor */
        ToolTest(const std::string&, const std::string&, const IInterface*);
        
        /** destructor */
        virtual ~ToolTest();
        
        /** AlgTool initilaize method */
        virtual StatusCode initialize() override;
        
        /** AlgTool finalize method*/
        virtual StatusCode finalize() override;
        
        virtual StatusCode toolTest() const override;

    private:
    
    ToolHandle<IToolTest2>    m_toolTest2;
                

};


#endif // ATS_GEOMETRYTOOLS_TOOLTEST_H

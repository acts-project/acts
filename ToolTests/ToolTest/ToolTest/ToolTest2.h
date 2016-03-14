///////////////////////////////////////////////////////////////////
// ToolTest2.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_GEOMETRYTOOLS_TOOLTEST2_H
#define ATS_GEOMETRYTOOLS_TOOLTEST2_H 1

// Core module
#include "CoreInterfaces/AlgToolBase.h"
// Geoemtry module
#include "ToolTestInterfaces/IToolTest2.h"
// Gaudi & Athena
#include "GaudiKernel/ToolHandle.h"

    
class ToolTest2 : public Ats::AlgToolBase, virtual public IToolTest2 {
        
    public:
        /** constructor */
        ToolTest2(const std::string&, const std::string&, const IInterface*);
        
        /** destructor */
        virtual ~ToolTest2();
        
        /** AlgTool initilaize method */
        virtual StatusCode initialize() override;
        
        /** AlgTool finalize method*/
        virtual StatusCode finalize() override;
        
        virtual StatusCode toolTest2() const override;

    private:
        
                

};


#endif // ATS_GEOMETRYTOOLS_TOOLTEST2_H

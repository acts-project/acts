///////////////////////////////////////////////////////////////////
// ToolTest3.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_GEOMETRYTOOLS_TOOLTEST3_H
#define ATS_GEOMETRYTOOLS_TOOLTEST3_H 1

// Core module
#include "CoreInterfaces/AlgToolBase.h"
// Geoemtry module
#include "ToolTestInterfaces/IToolTest3.h"
// Gaudi & Athena
#include "GaudiKernel/ToolHandle.h"

    
class ToolTest3 : public Ats::AlgToolBase, virtual public IToolTest3 {
        
    public:
        /** constructor */
        ToolTest3(const std::string&, const std::string&, const IInterface*);
        
        /** destructor */
        virtual ~ToolTest3();
        
        /** AlgTool initilaize method */
        virtual StatusCode initialize() override;
        
        /** AlgTool finalize method*/
        virtual StatusCode finalize() override;
        
        virtual StatusCode toolTest3() const override;

    private:
        
                

};


#endif // ATS_GEOMETRYTOOLS_TOOLTEST3_H

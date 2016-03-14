//////////////////////////////////////////////////////////////////
// TooAlgorithm.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_EXAMPLES_TOOLALGORITHM_H
#define ATS_EXAMPLES_TOOLALGORITHM_H

// Core
#include "CoreInterfaces/AlgorithmBase.h"
// Gauid Include
#include "GaudiKernel/ToolHandle.h"

#include "ToolTestInterfaces/IToolTest.h"
// Trk includes


class ToolAlgorithm : public Ats::AlgorithmBase  {
     public:

       /** Standard Athena-Algorithm Constructor */
       ToolAlgorithm(const std::string& name, ISvcLocator* pSvcLocator);

       /* initialize the test, i.e. retrieve the TrackingGeometry Svc */
       StatusCode initialize() override;
       
       /* exectue - and set executed to true */
       StatusCode execute() override;

       /* finalize  */
       StatusCode finalize() override;
              
       
     private:
        
    ToolHandle<IToolTest>    m_toolTest;
    
   };


#endif

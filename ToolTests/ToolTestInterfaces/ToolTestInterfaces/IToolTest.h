///////////////////////////////////////////////////////////////////
// IToolTest.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_GEOMETRYINTERFACES_ITOOLTEST_H
#define ATS_GEOMETRYINTERFACES_ITOOLTEST_H 1

// Gaudi
#include "GaudiKernel/IAlgTool.h"


  /** Interface ID for IToolTests*/  
  static const InterfaceID IID_IToolTest("IToolTest", 1, 0);

  class IToolTest : virtual public IAlgTool {
    
    public:
      /**Virtual destructor*/
      virtual ~IToolTest(){}
//       DeclareInterfaceID(IToolTest, 1, 0);
      /** AlgTool and IAlgTool interface methods */
      static const InterfaceID& interfaceID() { return IID_IToolTest; }

      /** Name identification */
      virtual StatusCode toolTest() const = 0;
             
  };

#endif // ATS_GEOMETRYINTERFACES_ITOOLTEST_H

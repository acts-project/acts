///////////////////////////////////////////////////////////////////
// IToolTest2.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_GEOMETRYINTERFACES_ITOOLTEST2_H
#define ATS_GEOMETRYINTERFACES_ITOOLTEST2_H 1

// Gaudi
#include "GaudiKernel/IAlgTool.h"


  /** Interface ID for IToolTest2s*/  
  static const InterfaceID IID_IToolTest2("IToolTest2", 1, 0);

  class IToolTest2 : virtual public IAlgTool {
    
    public:
      /**Virtual destructor*/
      virtual ~IToolTest2(){}
//       DeclareInterfaceID(IToolTest2, 1, 0);
      /** AlgTool and IAlgTool interface methods */
      static const InterfaceID& interfaceID() { return IID_IToolTest2; }

      /** Name identification */
      virtual StatusCode toolTest2() const = 0;
             
  };

#endif // ATS_GEOMETRYINTERFACES_ITOOLTEST2_H

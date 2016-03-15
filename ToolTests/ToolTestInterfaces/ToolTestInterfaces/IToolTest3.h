///////////////////////////////////////////////////////////////////
// IToolTest3.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_GEOMETRYINTERFACES_ITOOLTEST3_H
#define ATS_GEOMETRYINTERFACES_ITOOLTEST3_H 1

// Gaudi
#include "GaudiKernel/IAlgTool.h"


  /** Interface ID for IToolTest3s*/  
  static const InterfaceID IID_IToolTest3("IToolTest3", 1, 0);

  class IToolTest3 : virtual public IAlgTool {
    
    public:
      /**Virtual destructor*/
      virtual ~IToolTest3(){}
//       DeclareInterfaceID(IToolTest3, 1, 0);
      /** AlgTool and IAlgTool interface methods */
      static const InterfaceID& interfaceID() { return IID_IToolTest3; }

      /** Name identification */
      virtual StatusCode toolTest3() const = 0;
             
  };

#endif // ATS_GEOMETRYINTERFACES_ITOOLTEST3_H

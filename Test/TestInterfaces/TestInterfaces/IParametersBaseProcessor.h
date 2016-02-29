///////////////////////////////////////////////////////////////////
// IParametersBaseProcessor.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_TESTINTERFACES_IPARAMETERSBASEPROCESSOR_H
#define ATS_TESTINTERFACES_IPARAMETERSBASEPROCESSOR_H 1

// EventData module
#include "ParametersBase/TrackParametersBase.h"
// Gaudi 
#include "GaudiKernel/IAlgTool.h"
//STL
#include <string>

namespace Ats {

  /** Interface ID for IParametersBaseProcessors*/  
  static const InterfaceID IID_IParametersBaseProcessor("IParametersBaseProcessor", 1, 0);
  
  /** @class IParametersBaseProcessor
  
       Interface class IParametersBaseProcessors
    
      @author Andreas.Salzburger@cern.ch
    */
  class IParametersBaseProcessor : virtual public IAlgTool {
    
    public:
      /**Virtual destructor*/
      virtual ~IParametersBaseProcessor(){}
      
      /** AlgTool and IAlgTool interface methods */
      static const InterfaceID& interfaceID() { return IID_IParametersBaseProcessor; }

      /** Processor Action to work on ParameterBase   */
      virtual StatusCode process(const TrackParametersBase& ) = 0;

      /** Processor Action to work on ParameterBase vector  */
      virtual StatusCode process(const std::vector<const TrackParametersBase*>& pBaseVector) = 0;

      /** Processor Action initializations */
      virtual StatusCode initProcessor() = 0;

  };

} // end of namespace

#endif // ATS_TESTINTERFACES_IPARAMETERSBASEPROCESSOR_H

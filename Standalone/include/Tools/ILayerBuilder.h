///////////////////////////////////////////////////////////////////
// ILayerBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYINTERFACES_ILAYERBUILDER_H
#define ACTS_GEOMETRYINTERFACES_ILAYERBUILDER_H 1

// Gaudi
#include "GaudiKernel/IAlgTool.h"
// STL
#include <vector>
#include <string>

namespace Acts {

  class Layer;
  typedef std::shared_ptr<const Layer> LayerPtr;
  typedef std::vector< LayerPtr > LayerVector;

  /** Interface ID for ILayerBuilders*/  
  static const InterfaceID IID_ILayerBuilder("ILayerBuilder", 1, 0);
  
  /** @class ILayerBuilder
  
    Interface class for ILayerBuilders in a typical 
    | EC- | Central | EC+ | 
    detector setup.
      
    @author Andreas.Salzburger@cern.ch
    */
  class ILayerBuilder : virtual public IAlgTool {
    
    public:
      /**Virtual destructor*/
      virtual ~ILayerBuilder(){}
//       DeclareInterfaceID(ILayerBuilder, 1, 0);
      /** AlgTool and IAlgTool interface methods */
      static const InterfaceID& interfaceID() { return IID_ILayerBuilder; }

      /** LayerBuilder interface method - returning the layers at negative side */
      virtual const LayerVector* negativeLayers() const = 0; 
      
      /** LayerBuilder interface method - returning the central layers */
      virtual const LayerVector* centralLayers() const = 0; 
      
      /** LayerBuilder interface method - returning the layers at negative side */
      virtual const LayerVector* positiveLayers() const = 0; 

      /** Name identification */
      virtual const std::string& identification() const = 0;
             
  };

} // end of namespace

#endif // ACTS_GEOMETRYINTERFACES_ILAYERBUILDER_H

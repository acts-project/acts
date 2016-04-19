///////////////////////////////////////////////////////////////////
// ILayerBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYINTERFACES_ILAYERBUILDER_H
#define ACTS_GEOMETRYINTERFACES_ILAYERBUILDER_H 1

// STL
#include <vector>
#include <string>
#include <memory>

namespace Acts {

  class Layer;
  typedef std::shared_ptr<const Layer> LayerPtr;
  typedef std::vector< LayerPtr > LayerVector;

  /** @class ILayerBuilder
  
    Interface class for ILayerBuilders in a typical 
    | EC- | Central | EC+ | 
    detector setup.
      
    @author Andreas.Salzburger@cern.ch
    */
  class ILayerBuilder
  { 
    public:
      /**Virtual destructor*/
      virtual ~ILayerBuilder(){}

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

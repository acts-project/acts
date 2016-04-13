///////////////////////////////////////////////////////////////////
// ILayerArrayCreator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYINTERFACES_ILAYERARRAYCREATOR_H
#define ACTS_GEOMETRYINTERFACES_ILAYERARRAYCREATOR_H 1

// Gaudi
#include "GaudiKernel/IAlgTool.h"
// Geometry module
#include "GeometryUtils/BinnedArray.h"
#include "GeometryUtils/BinningType.h"
// STL
#include <vector>

namespace Acts {

  /** forward declarations*/
  class Layer;
  typedef std::shared_ptr<const Layer> LayerPtr;

  /** @typedef LayerArray and LayerVector */
  typedef BinnedArray< LayerPtr > LayerArray;
  typedef std::vector< LayerPtr > LayerVector;
  
  /** Interface ID for ILayerArrayCreators*/  
  static const InterfaceID IID_ILayerArrayCreator("ILayerArrayCreator", 1, 0);
  
  /** @class ILayerArrayCreator
  
    Interface class ILayerArrayCreators, it inherits from IAlgTool. 
    
    It receives the LayerVector and creaets an array with NaivgationLayer objects
    filled in between.
    
    @author Andreas.Salzburger@cern.ch
  */
  
  class ILayerArrayCreator : virtual public IAlgTool {
    
    public:
      /**Virtual destructor*/
      virtual ~ILayerArrayCreator(){}
      
    //   DeclareInterfaceID(ILayerArrayCreator, 1, 0);
      /** AlgTool and IAlgTool interface methods */
      static const InterfaceID& interfaceID() { return IID_ILayerArrayCreator; }

      /** LayerArraycreator interface method */
      virtual LayerArray* layerArray(const LayerVector& layers, 
                                     double min,
                                     double max,
                                     BinningType btype = arbitrary,
                                     BinningValue bv   = binX) const = 0; 
   
  };


} // end of namespace


#endif // ACTS_GEOMETRYINTERFACES_ILAYERARRAYCREATOR_H



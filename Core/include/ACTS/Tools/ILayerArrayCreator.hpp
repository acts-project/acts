///////////////////////////////////////////////////////////////////
// ILayerArrayCreator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYINTERFACES_ILAYERARRAYCREATOR_H
#define ACTS_GEOMETRYINTERFACES_ILAYERARRAYCREATOR_H 1

// Geometry module
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Utilities/BinningType.hpp"
// STL
#include <vector>
#include <memory>

namespace Acts
{
  /** forward declarations*/
  class Layer;
  typedef std::shared_ptr<const Layer> LayerPtr;

  /** @typedef LayerArray and LayerVector */
  typedef BinnedArray< LayerPtr > LayerArray;
  typedef std::vector< LayerPtr > LayerVector;
  
  /** @class ILayerArrayCreator
  
    Interface class ILayerArrayCreators, it inherits from IAlgTool. 
    
    It receives the LayerVector and creaets an array with NaivgationLayer objects
    filled in between.
    
    @author Andreas.Salzburger@cern.ch
  */
  
  class ILayerArrayCreator
  {  
    public:
      /**Virtual destructor*/
    virtual ~ILayerArrayCreator() = default;
      
      /** LayerArraycreator interface method */
      virtual std::unique_ptr<const LayerArray> layerArray(const LayerVector& layers,
                                                           double min,
                                                           double max,
                                                           BinningType btype = arbitrary,
                                                           BinningValue bv   = binX) const = 0; 
   
  };
} // end of namespace


#endif // ACTS_GEOMETRYINTERFACES_ILAYERARRAYCREATOR_H

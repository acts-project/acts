///////////////////////////////////////////////////////////////////
// LayerArrayCreator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYTOOLS_LAYERARRAYCREATOR_H
#define ACTS_GEOMETRYTOOLS_LAYERARRAYCREATOR_H 1

#ifndef TRKDETDESCR_TAKESMALLERBIGGER
#define TRKDETDESCR_TAKESMALLERBIGGER
#define takeSmaller(current,test) current = current < test ? current : test
#define takeBigger(current,test)  current = current > test ? current : test
#define takeSmallerBigger(cSmallest, cBiggest, test) takeSmaller(cSmallest, test); takeBigger(cBiggest, test)
#endif

// Core module
#include "ACTS/Utilities/Definitions.h"
// Geometry module
#include "ACTS/Tools/ILayerArrayCreator.h"
// STL
#include <algorithm>

namespace Acts {

    class Surface;
    class Layer;

    /** @class LayerArrayCreator

      The LayerArrayCreator is a simple Tool that helps to construct
      LayerArrays from std::vector of Acts::CylinderLayer, Acts::DiscLayer, Acts::PlaneLayer.

      It fills the gaps automatically with Acts::NavigationLayer to be processed easily in the
      Navigation of the Extrapolation process.

     @TODO Julia: make private tools private again after Gaudi update (bug in Gaudi), marked with //b
     
      @author Andreas.Salzburger@cern.ch   
     */

    class LayerArrayCreator : public ILayerArrayCreator {

      public:
        /** Constructor */
        LayerArrayCreator() = default;
        
        /** Destructor */
        virtual ~LayerArrayCreator() = default;

        /** LayerArraycreator interface method 
           - we assume the layer thickness to be used together with the binning value */
        LayerArray* layerArray(const LayerVector& layers, 
                               double min,
                               double max,
                               BinningType btype = arbitrary,
                               BinningValue bvalue = binX) const override; 
      
      private:
          Surface* createNavigationSurface(const Layer& layer, BinningValue bvalue, double offset) const;
    };

} // end of namespace

#endif // ACTS_GEOMETRYTOOLS_LAYERARRAYCREATOR_H


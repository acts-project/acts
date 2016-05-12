///////////////////////////////////////////////////////////////////
// LayerArrayCreator.cpp, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "ACTS/Tools/LayerArrayCreator.hpp"

#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/BinnedArray1D.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Utilities/GeometryStatics.hpp"
#include "ACTS/Utilities/GeometryObjectSorter.hpp"
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Layers/NavigationLayer.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Surfaces/DiscSurface.hpp"
#include "ACTS/Surfaces/DiscBounds.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/SurfaceBounds.hpp"
#include "ACTS/Surfaces/TrapezoidBounds.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"

Acts::LayerArray* Acts::LayerArrayCreator::layerArray(const LayerVector& layersInput,
                                                      double min, double max, 
                                                      BinningType bType,
                                                      BinningValue bValue) const 
{

   // MSG_VERBOSE( "Build LayerArray with "     << layersInput.size() << " layers at input." );
   // MSG_VERBOSE( "       min/max provided : " << min << " / " << max );
   // MSG_VERBOSE( "       binning type     : " << bType );
   // MSG_VERBOSE( "       binning value    : " << bValue );

   // create a local copy of the layer vector
   LayerVector layers(layersInput);

   // sort it accordingly to the binning value
   GeometryObjectSorterT< std::shared_ptr<const Layer> > layerSorter(bValue);
   std::sort(layers.begin(),layers.end(),layerSorter);
   // useful typedef
   typedef std::pair< std::shared_ptr<const Layer>, Vector3D > LayerOrderPosition;
   // needed for all cases
   std::shared_ptr<const Layer>   layer      = nullptr;
   LayerArray*                    layerArray = nullptr;
   BinUtility*                    binUtility = nullptr;
   std::vector< LayerOrderPosition >   layerOrderVector;

   // switch the binning type
   switch (bType) {

        // equidistant binning - no navigation layers built - only equdistant layers
        case equidistant :
        {
            // loop over layers and put them in
            for (auto& layIter : layers ) {
                // MSG_VERBOSE( "equidistant : registering a Layer at binning position : " << toString(layIter->binningPosition(bValue)) );
                layerOrderVector.push_back( LayerOrderPosition(layIter, layIter->binningPosition(bValue) ));
            }
            // create the binUitlity
            binUtility = new BinUtility(layers.size(), min, max, open, bValue);
            // MSG_VERBOSE( "equidistant : created a BinUtility as " << *binUtility );
        } break;

        // arbitrary binning
        case arbitrary :
        {
            std::vector<float> boundaries;
            // initial step
            boundaries.push_back(min);
            double layerValue      = 0.;
            double layerThickness  = 0.;
            std::shared_ptr<const Layer> navLayer  = nullptr;
            std::shared_ptr<const Layer> lastLayer = nullptr;
            // loop over layers
            for (auto& layIter : layers ) {
                // estimate the offset
                layerThickness  = layIter->thickness();
                layerValue      = layIter->binningPositionValue(bValue);
                // register the new boundaries in the step vector
                boundaries.push_back(layerValue-0.5*layerThickness);
                boundaries.push_back(layerValue+0.5*layerThickness);
                // calculate the layer value for the offset
                double navigationValue = 0.5*((layerValue-0.5*layerThickness) + boundaries[boundaries.size()-3]);
                //if layers are attached to each other, no navigation layer needed
                if (navigationValue!=(layerValue-0.5*layerThickness)) {
                    // create the navigation layer surface from the layer
                    Surface* navLayerSurface = createNavigationSurface(*layIter, bValue, -fabs(layerValue-navigationValue));
                    navLayer = NavigationLayer::create(navLayerSurface);
                    // push the navigation layer in
                    layerOrderVector.push_back(LayerOrderPosition(navLayer, navLayer->binningPosition(bValue)));
                    // MSG_VERBOSE( "arbitrary : creating a  NavigationLayer at " << toString(navLayerSurface->binningPosition(bValue)) );
                }
                // push the original layer in
                layerOrderVector.push_back(LayerOrderPosition(layIter, layIter->binningPosition(bValue) ));
                // MSG_VERBOSE( "arbitrary : registering MaterialLayer at  " <<  toString(layIter->binningPosition(bValue)));
                // remember the last
                lastLayer = layIter;
            }
            // a final navigation layer
            // calculate the layer value for the offset
            double navigationValue = 0.5*(boundaries[boundaries.size()-1]+max);
            //create navigation layer only when necessary
            if (navigationValue!=max) {
                // create the navigation layer surface from the layer
                Surface* navLayerSurface = createNavigationSurface(*lastLayer, bValue, navigationValue-layerValue);
                navLayer = NavigationLayer::create(navLayerSurface);
                // push the navigation layer in
                layerOrderVector.push_back(LayerOrderPosition(navLayer, navLayer->binningPosition(bValue)));
                // MSG_VERBOSE( "arbitrary : creating a  NavigationLayer at " << toString(navLayerSurface->binningPosition(bValue)) );
            }
            // now close the boundaries
            boundaries.push_back(max);
            // some screen output
            // MSG_VERBOSE( layerOrderVector.size() << " Layers (material + navigation) built. " );
            // create the BinUtility
            binUtility = new BinUtility(boundaries, open, bValue);
            // MSG_VERBOSE( "arbitrary : created a BinUtility as " << *binUtility );

        } break;
        // default return nullptr
        default : { return nullptr; }

    }
    // create the BinnedArray
    layerArray = new BinnedArray1D< std::shared_ptr< const Layer> >(layerOrderVector, binUtility);
    // return what we have here
    return layerArray;

}

Acts::Surface* Acts::LayerArrayCreator::createNavigationSurface(const Layer& layer, BinningValue bValue, double offset) const
{
    // surface reference
    const Surface& layerSurface = layer.surfaceRepresentation();
    // translation to be applied
    Vector3D translation(0.,0.,0.);
    // switching he binnig values
    switch (bValue) {
        // case x
        case binX : { translation = Vector3D(offset, 0., 0.);} break;
        // case y
        case binY : { translation = Vector3D(0., offset, 0.);} break;
        // case z
        case binZ : { translation = Vector3D(0., 0., offset); } break;
        // case R
        case binR : {
            // binning in R and cylinder surface means something different
            if ( layerSurface.type() == Surface::Cylinder ) break;
            translation = Vector3D(offset, 0., 0.);
        } break;
        // do nothing for the default
        default : {
            // MSG_WARNING("Not yet implemented.");
        }
    }
    // navigation surface
    Surface* navigationSurface = nullptr;
    // for everything else than a cylinder it's a copy with shift
    if (layerSurface.type() != Surface::Cylinder) {
        // create a transform that does the shift
        Transform3D* shift = new Transform3D;
        (*shift) = Translation3D(translation);
        navigationSurface = layerSurface.clone(shift);
        // delete the shift again
        delete shift;
    } else {
       // get the bounds
       const CylinderBounds* cBounds = dynamic_cast<const CylinderBounds*>(&(layerSurface.bounds()));
       double navigationR = cBounds->r()+offset;
       double halflengthZ = cBounds->halflengthZ();
       // new navigation layer
       navigationSurface = new CylinderSurface(layerSurface.cachedTransform(), navigationR, halflengthZ);
    }
    return navigationSurface;
}

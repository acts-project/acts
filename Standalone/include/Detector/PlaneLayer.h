///////////////////////////////////////////////////////////////////
// PlaneLayer.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DETECTOR_PLANELAYER_H
#define ACTS_DETECTOR_PLANELAYER_H

class MsgStream;

// Geomety module
#include "Detector/Layer.h"
#include "Surfaces/PlaneSurface.h"
// Core module
#include "EventDataUtils/PropDirection.h"
// STL sorting
#include <algorithm>

namespace Acts {

  class OverlapDescriptor;
  class ApproachDescriptor;

  /** 
   @class PlaneLayer
   
   Class to describe a planar detector layer for tracking, 
   it inhertis from both, Layer base class and PlaneSurface class
       
   @author Andreas.Salzburger@cern.ch 
   */

  class PlaneLayer : virtual public PlaneSurface, public Layer {
      
      public:
        /** Factory for a shared plane layer */
        static LayerPtr create(std::shared_ptr<Transform3D> transform,
                                                   std::shared_ptr<const PlanarBounds> pbounds,
                                                   SurfaceArray* surfaces = nullptr,
                                                   double thickness = 0.,
                                                   OverlapDescriptor* od = nullptr,
                                                   ApproachDescriptor* ad = nullptr,
                                                   int laytyp=int(Acts::active)) 
        { return LayerPtr(new PlaneLayer(transform, pbounds, surfaces, thickness, od, ad, laytyp)); } 
                           
        /** Factory for a shared plane layer */
        static LayerPtr create(const PlaneLayer& pla, const Transform3D& tr)
        { return LayerPtr(new PlaneLayer(pla,tr)); }  
        
        /** Clone with a shift - only cloning that is allowed */
        LayerPtr cloneWithShift(const Transform3D& shift) const 
        { return PlaneLayer::create(*this,shift); }                 
                           
        /** Copy constructor of PlaneLayer - forbidden */
        PlaneLayer(const PlaneLayer& pla) = delete;
        
        /** Assignment operator for PlaneLayers - forbidden */
        PlaneLayer& operator=(const PlaneLayer&) = delete;
               
        /** Destructor*/
        virtual ~PlaneLayer(){}  
    
        /** Transforms the layer into a Surface representation for extrapolation */
        const PlaneSurface& surfaceRepresentation() const;            
    
     private:

        /** build approach surfaces */
        void buildApproachDescriptor() const;

    protected:
        /** Default Constructor*/
        PlaneLayer(){}
                   
        /** Constructor with PlaneSurface components  
           - shared bounds */
        PlaneLayer(std::shared_ptr<Transform3D> transform,
                   std::shared_ptr<const PlanarBounds>& pbounds,
                   SurfaceArray* = nullptr,
                   double thickness = 0.,
                   OverlapDescriptor* od = nullptr,
                   ApproachDescriptor* ad = nullptr,
                   int laytyp=int(Acts::active));
                   
        /** Copy constructor with shift*/
        PlaneLayer(const PlaneLayer& pla, const Transform3D& tr);

  };

} // end of namespace

#endif // TRKGEOMETY_PLANELAYER_H

///////////////////////////////////////////////////////////////////
// SubtractedPlaneLayer.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DETECTOR_SUBTRACTEDPLANELAYER_H
#define ACTS_DETECTOR_SUBTRACTEDPLANELAYER_H

// Geometry module
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Layers/PlaneLayer.hpp"
#include "ACTS/Surfaces/SubtractedPlaneSurface.hpp"
// Core module

namespace Acts {

  class OverlapDescriptor;

  /** 
   @class SubtractedPlaneLayer
   
   Class to describe a planar detector layer for tracking, with subtraction 
       
   @author Sarka.Todorova@cern.ch, Andres.Salzburger@cern.ch 
   */

  class SubtractedPlaneLayer : virtual public SubtractedPlaneSurface, public Layer {
      
      public:
        /** Factory Constructor with SubtractedPlaneSurface  */
        static LayerPtr create(const SubtractedPlaneSurface* subtrPlaneSurf,
                                                   double thickness = 0.,
                                                   OverlapDescriptor* od = 0,
                                                   int laytyp=int(Acts::active))
        { return LayerPtr(new SubtractedPlaneLayer(subtrPlaneSurf,thickness, od, laytyp)); }
        
        /** Factory Constructor as a copy with shift */
        static LayerPtr create(const SubtractedPlaneLayer& pla, const Transform3D& tr)
        { return LayerPtr(new SubtractedPlaneLayer(pla,tr)); }
            
        /** Clone with shift - the only allowed way to clone */
        LayerPtr cloneWithShift(const Transform3D& shift) const override
        { return SubtractedPlaneLayer::create(*this,shift); }      
                           
        /** Copy constructor of SubtractedPlaneLayer - is forbidden */
        SubtractedPlaneLayer(const SubtractedPlaneLayer& pla) = delete;
        
        /** Assignment operator for PlaneLayers */
        SubtractedPlaneLayer& operator=(const SubtractedPlaneLayer&) = delete;
                                          
        /** Destructor*/
        virtual ~SubtractedPlaneLayer(){}   
    
        /** Transforms the layer into a Surface representation for extrapolation */
        const SubtractedPlaneSurface& surfaceRepresentation() const override;            
        
    protected:
        /** Default Constructor*/
        SubtractedPlaneLayer(){}
        
        /** Constructor with SubtractedPlaneSurface and Material  */
        SubtractedPlaneLayer(const SubtractedPlaneSurface* subtrPlaneSurf,
                             double thickness = 0.,
                             OverlapDescriptor* od = 0,
                             int laytyp=int(Acts::active)); 
                             
        /** Copy constructor with shift*/
        SubtractedPlaneLayer(const SubtractedPlaneLayer& pla, const Transform3D& tr);                             
                                                                   
  };

} // end of namespace

#endif // TRKGEOMETY_SUBTRACTEDPLANELAYER_H


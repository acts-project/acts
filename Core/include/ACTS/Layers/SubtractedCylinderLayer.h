///////////////////////////////////////////////////////////////////
// SubtractedCylinderLayer.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DETECTOR_SUBTRACTEDCYLINDERLAYER_H
#define ACTS_DETECTOR_SUBTRACTEDCYLINDERLAYER_H

// Trk
#include "ACTS/Utilities/Definitions.h"
#include "ACTS/Layers/Layer.h"
#include "ACTS/Surfaces/SubtractedCylinderSurface.h"
// Amg

namespace Acts {

  /**
   @class SubtractedCylinderLayer
   
   Class to describe a cylindrical detector layer for tracking, with subtraction; it inhertis from both, 
   Layer base class and SubtractedCylinderSurface class
       
   @author Sarka.Todorova@cern.ch
  */

  class SubtractedCylinderLayer : virtual public SubtractedCylinderSurface, public Layer {
                   
      public:
        /** Factory constructor with arguments */
        static LayerPtr create(const SubtractedCylinderSurface* subCyl,
                                                   double thickness = 0.,
                                                   int laytyp=int(Acts::active))
        { return LayerPtr(new SubtractedCylinderLayer(subCyl,thickness, laytyp));} 
          
        /** Factory copy constructor with shift */  
        static LayerPtr create(const SubtractedCylinderLayer& cla, const Transform3D& tr)
        { return LayerPtr(new SubtractedCylinderLayer(cla,tr));}          
                              
        /** Clone with shift - the only allowed way to clone */
        LayerPtr cloneWithShift(const Transform3D& shift) const override
        { return SubtractedCylinderLayer::create(*this,shift); }                   
                              
        /** Copy constructor - forbidden */
        SubtractedCylinderLayer(const SubtractedCylinderLayer&) = delete;
        
        /** Assignment operator - forbidden */
        SubtractedCylinderLayer& operator=(const SubtractedCylinderLayer&) = delete;
                      
        /** Destructor*/
        virtual ~SubtractedCylinderLayer(){}  
        
        /** Transforms the layer into a Surface representation for extrapolation */
        const SubtractedCylinderSurface& surfaceRepresentation() const override;
        
        /** use the base class insideBounds (Vector2d, BoundaryCheck) */
        using CylinderSurface::insideBounds;


   private:
       /** Default Constructor*/
       SubtractedCylinderLayer(){}
       
       /** Constructor with SubtractedCylinderSurface components */
       SubtractedCylinderLayer(const SubtractedCylinderSurface* subCyl,
                               double thickness = 0.,
                               int laytyp=int(Acts::active));       

       /** Copy constructor with shift*/
       SubtractedCylinderLayer(const SubtractedCylinderLayer& cla, const Transform3D& tr);
 
  };
 
} // end of namespace

#endif // ACTS_DETECTOR_SUBTRACTEDCYLINDERLAYER_H


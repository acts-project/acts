///////////////////////////////////////////////////////////////////
// ConeLayer.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DETECTOR_CONELAYER_H
#define ACTS_DETECTOR_CONELAYER_H

class MsgStream;

// Core module
#include "ACTS/Utilities/Definitions.hpp"
// Geometry module
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Surfaces/ConeSurface.hpp"
// STL sorting
#include <algorithm>

namespace Acts {
  
  class ConeBounds;
  class OverlapDescriptor;
  
  /**
     @class ConeLayer
     
     Class to describe a cylindrical detector layer for tracking, it inhertis from both, 
     Layer base class and ConeSurface class
     
     @author Ian.Watson@cern.ch, Andreas.Salzburger@cern.ch
  */
  
  class ConeLayer : virtual public ConeSurface, public Layer {
    
    public:
      /** Factory for shared layer */
      static LayerPtr create(std::shared_ptr<Transform3D> transform,
                                                 std::shared_ptr<const ConeBounds> cbounds,
                                                 SurfaceArray* surfaceArray,
                                                 double thickness = 0.,
                                                 OverlapDescriptor* od = 0,
                                                 int laytyp=int(Acts::active)) 
      { return LayerPtr(new ConeLayer(transform,cbounds,surfaceArray, thickness, od, laytyp));}

      /** Factory for shared layer with shift */
      static LayerPtr create(const ConeLayer& cla, const Transform3D& shift)
      { return LayerPtr(new ConeLayer(cla, shift)); }  
      
      /** Clone with a shift - only cloning that is allowed */
      LayerPtr cloneWithShift(const Transform3D& shift) const override
      { return ConeLayer::create(*this,shift); }
      
      /** Copy constructor of ConeLayer - forbidden */
      ConeLayer(const ConeLayer& cla) = delete;
    
      /** Assignment operator for ConeLayers - forbidden */
      ConeLayer& operator=(const ConeLayer&) = delete;
    
      /** Destructor*/
      virtual ~ConeLayer(){}  
    
      /** Transforms the layer into a Surface representation for extrapolation */
      const ConeSurface& surfaceRepresentation() const override;
            
    protected:  
      /** Default Constructor*/
      ConeLayer(){}
      
      /** Constructor with ConeSurface components,
         MaterialProperties and pointer SurfaceArray (passing ownership) 
         @TODO implement ApproachDescriptor */
      ConeLayer(std::shared_ptr<Transform3D> transform,
                std::shared_ptr<const ConeBounds> cbounds,
                SurfaceArray* surfaceArray,
                double thickness = 0.,
                OverlapDescriptor* od = 0,
                int laytyp=int(Acts::active));
                
       /** Copy constructor with shift*/
       ConeLayer(const ConeLayer& cla, const Transform3D& tr);                
         
  };

} // end of namespace

#endif // TRKGEOMETY_CONELAYER_H

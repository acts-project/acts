///////////////////////////////////////////////////////////////////
// DiscLayer.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DETECTOR_DISCLAYER_H
#define ACTS_DETECTOR_DISCLAYER_H 1

class MsgStream;

#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Surfaces/DiscSurface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
// STL sorting
#include <algorithm>

namespace Acts {

  class DiscBounds;
  class SurfaceMaterial;
  class OverlapDescriptor;
  class ApproachDescriptor;
  
  /**
   @class DiscLayer
   
   Class to describe a disc-like detector layer for tracking, 
   it inhertis from both, Layer base class
   and DiscSurface class
       
   @author Andreas.Salzburger@cern.ch 
   */

  class DiscLayer : virtual public DiscSurface, public Layer {

     friend class TrackingVolume;
      
     public:                  
       /** Factory constructor with DiscSurface components and pointer to SurfaceArray (passing ownership) */
       static LayerPtr create(std::shared_ptr<Transform3D> transform,
                              std::shared_ptr<const DiscBounds> dbounds,
                              std::unique_ptr<SurfaceArray> surfaceArray = nullptr,
                              double thickness = 0.,
                              OverlapDescriptor* od = nullptr,
                              ApproachDescriptor* ad = nullptr,
                              int laytyp=int(Acts::passive)) 
      { return LayerPtr(new DiscLayer(transform,dbounds,std::move(surfaceArray),thickness,od,ad,laytyp)); }
                                                        
       /** Factory constructor as copy with shift */
       static LayerPtr create(const DiscLayer& cla, const Transform3D& shift)
       { return LayerPtr(new DiscLayer(cla,shift)); }                                                 

       /** Clone with a shift - only cloning that is allowed */
       LayerPtr cloneWithShift(const Transform3D& shift) const override
       { return DiscLayer::create(*this,shift); }

       /** Copy constructor of DiscLayer - forbidden */
       DiscLayer(const DiscLayer& cla) = delete;
       
       /* *Assignment operator for DiscLayers - forbidden */
       DiscLayer& operator=(const DiscLayer&) = delete;
             
       /** Destructor*/
       virtual ~DiscLayer(){}  
               
       /** Transforms the layer into a Surface representation for extrapolation */
       const DiscSurface& surfaceRepresentation() const override;
     
     private:   
       /** build approach surfaces */
       void buildApproachDescriptor() const;
    
     protected:
       /** Default Constructor*/
       DiscLayer(){}

       /** Constructor with DiscSurface components and pointer to SurfaceArray (passing ownership) */
       DiscLayer(std::shared_ptr<Transform3D> transform,
                 std::shared_ptr<const DiscBounds> dbounds,
                 std::unique_ptr<SurfaceArray> surfaceArray = nullptr,
                 double thickness = 0.,
                 OverlapDescriptor* od = nullptr,
                 ApproachDescriptor* ad = nullptr,
                 int laytyp=int(Acts::active));   
                 
       /** Copy constructor with shift*/
       DiscLayer(const DiscLayer& cla, const Transform3D& tr);              
  };

} // end of namespace

#endif // TRKGEOMETY_DISCLAYER_H

///////////////////////////////////////////////////////////////////
// ITrackingVolumeBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYINTERFACES_ITRACKINGVOLUMEBUILDER_H
#define ACTS_GEOMETRYINTERFACES_ITRACKINGVOLUMEBUILDER_H 1

// Gaudi
#include "GaudiKernel/IAlgTool.h"
#include <tuple>

namespace Acts {

  class VolumeBounds;
  class TrackingVolume;
  class Layer;
  class Volume;
  typedef std::shared_ptr<const TrackingVolume> TrackingVolumePtr;
  typedef std::shared_ptr<const VolumeBounds>   VolumeBoundsPtr;
  typedef std::shared_ptr<const Volume>         VolumePtr;
  typedef std::tuple< std::vector< std::shared_ptr<const Layer> >*, std::vector< std::shared_ptr<const Layer> >*, std::vector< std::shared_ptr<const Layer> >* > LayerTriple;
    typedef std::tuple< VolumePtr, VolumePtr, VolumePtr> VolumeTriple;
      
  
  /** Interface ID for ITrackingVolumeBuilders*/  
  static const InterfaceID IID_ITrackingVolumeBuilder("ITrackingVolumeBuilder", 1, 0);
  
  /** @class ITrackingVolumeBuilder
  
     Interface class ITrackingVolumeBuilders & inherits from IAlgTool

     this returns the sub-detector tracking volume that is wrapped by the next outer one
     in the TrackingGeometry building process
  
     If an innerVolume is given, this is wrapped
     If a VolumeBounds object is given this defines the maximum extent.
  
     @author Andreas.Salzburger@cern.ch
    */
  class ITrackingVolumeBuilder : virtual public IAlgTool {
    
    public:
      /**Virtual destructor*/
      virtual ~ITrackingVolumeBuilder(){}
      
//       DeclareInterfaceID(ITrackingVolumeBuilder, 1, 0);
      /** AlgTool and IAlgTool interface methods */
      static const InterfaceID& interfaceID() { return IID_ITrackingVolumeBuilder; }

      /** TrackingVolumeBuilder interface method - returns the volumes of Volumes */
      virtual TrackingVolumePtr trackingVolume(TrackingVolumePtr insideVolume = nullptr,
                                               VolumeBoundsPtr outsideBounds = nullptr,
                                               const LayerTriple* layerTriple = nullptr,
                                               const VolumeTriple* volumeTriple = nullptr) const = 0;
      
                       
  };

} // end of namespace

#endif // ACTS_GEOMETRYINTERFACES_IITRACKINGVOLUMEBUILDER_H

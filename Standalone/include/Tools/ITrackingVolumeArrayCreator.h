///////////////////////////////////////////////////////////////////
// ITrackingVolumeArrayCreator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYINTERFACES_ITRACKINGVOLUMEARRAYCREATOR_H
#define ACTS_GEOMETRYINTERFACES_ITRACKINGVOLUMEARRAYCREATOR_H 1

// Gaudi
#include "GaudiKernel/IAlgTool.h"
// Geometry module
#include "GeometryUtils/BinnedArray.h"
#include "GeometryUtils/BinningType.h"
// STL
#include <vector>

namespace Acts {

  /** forward declarations*/
  class TrackingVolume;
  typedef std::shared_ptr< const TrackingVolume > TrackingVolumePtr;

  /** @typedef TrackingVolumeArray */
  typedef BinnedArray< TrackingVolumePtr > TrackingVolumeArray;
  typedef std::vector< TrackingVolumePtr > TrackingVolumeVector;
  
  /** Interface ID for ITrackingVolumeArrayCreator*/  
  static const InterfaceID IID_ITrackingVolumeArrayCreator("ITrackingVolumeArrayCreator", 1, 0);
  
  /** @class ITrackingVolumeArrayCreator
    
    Interface class ITrackingVolumeArrayCreators It inherits from IAlgTool. 
    
    It is designed to centralize the code to create
    Arrays of Tracking Volumes for both:
  
      - confinement in another TrackingVolume
      - navigation and glueing
  
    Arrays for glueing and confinement are often the same, 
    therefore the newly created TrackingVolumeArray is done by a shared_ptr

    @author Andreas.Salzburger@cern.ch
    */
  class ITrackingVolumeArrayCreator : virtual public IAlgTool {
    
    public:
      /**Virtual destructor*/
      virtual ~ITrackingVolumeArrayCreator(){}
      
//       DeclareInterfaceID(ITrackingVolumeArrayCreator, 1, 0);
      /** AlgTool and IAlgTool interface methods */
      static const InterfaceID& interfaceID() { return IID_ITrackingVolumeArrayCreator; }

      /** TrackingVolumeArrayCreator interface method - creates array depending on the binning type */
      virtual std::shared_ptr<const TrackingVolumeArray> trackingVolumeArray(const TrackingVolumeVector& vols, BinningValue bVal) const = 0; 
      
  };

} // end of namespace

#endif


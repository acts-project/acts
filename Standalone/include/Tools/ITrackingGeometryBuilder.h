///////////////////////////////////////////////////////////////////
// ITrackingGeometryBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYINTERFACES_ITRACKINGGEOMETRYBUILDER_H
#define ACTS_GEOMETRYINTERFACES_ITRACKINGGEOMETRYBUILDER_H 1

// Gaudi
#include "GaudiKernel/IAlgTool.h"

namespace Acts {

  class TrackingGeometry;

  /** Interface ID for ITrackingGeometryBuilders*/  
  static const InterfaceID IID_ITrackingGeometryBuilder("ITrackingGeometryBuilder", 1, 0);
   // DeclareInterfaceID(ITrackingGeometryBuilder, 1, 0);
  
  /** @class ITrackingGeometryBuilder
    
    Interface class for the TrackingGeometry building,
    this is used by the TrackingGeometrySvc to build the geoemtry.
  
    The TrackingGeometry is written to the detector store and thus not created
    as a std::shared_ptr.
  
    The TrackingGeometry is returned as a non-const object in order to recreate
    from conditions callback if necessary.
      
    @author Andreas.Salzburger@cern.ch
    */
  class ITrackingGeometryBuilder : virtual public IAlgTool {
    
    public:
//      DeclareInterfaceID(ITrackingGeometryBuilder, 1, 0);
      /**Virtual destructor*/
      virtual ~ITrackingGeometryBuilder(){}
      
      /** AlgTool and IAlgTool interface methods */
      static const InterfaceID& interfaceID() { return IID_ITrackingGeometryBuilder; }

      /** TrackingGeometry Interface methode */
      virtual TrackingGeometry* trackingGeometry() const = 0;
      
  };

} // end of namespace

#endif // ACTS_GEOMETRYINTERFACES_IGEOMETRYBUILDER_H

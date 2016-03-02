///////////////////////////////////////////////////////////////////
// IDD4hepLayerHelper.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_DD4HEPGEOMETRYINTERFACES_IDD4HEPLAYERHELPER_H
#define ATS_DD4HEPGEOMETRYINTERFACES_IDD4HEPLAYERHELPER_H 1

// Gaudi
#include "GaudiKernel/IAlgTool.h"
#include "GeometryInterfaces/ITrackingVolumeBuilder.h"

namespace DD4hep {
    namespace Geometry {
        class DetElement;
    }
}

namespace Add4hep {
    
    static const InterfaceID IID_IDD4hepLayerHelper("IDD4hepLayerHelper", 1, 0);
    
    /** @ class IDD4hepLayerHelper
     
     Interface to create a triplet of layers out of a DD4hep::DetElement.
     
     @author julia.hrdinka@aon.at
     */
    
    class IDD4hepLayerHelper : virtual public IAlgTool {
    
    public:
        /** Virtual destructor */
        virtual ~IDD4hepLayerHelper(){}
        
        /** AlgTool and IAlgTool interface methods */
        static const InterfaceID& interfaceID() { return IID_IDD4hepLayerHelper; }
        
        /** Creates a triple of all the three possible Layer types of the given volume detector element*/
        virtual const Ats::LayerTriple* createLayerTriple(DD4hep::Geometry::DetElement& motherDetElement) = 0;
    };
}

#endif //ATS_DD4HEPGEOMETRYINTERFACES_IDD4HEPLAYERHELPER_H

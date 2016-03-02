///////////////////////////////////////////////////////////////////
// LayerArrayCreator.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_DD4HEPGEOMETRYINTERFACES_IDD4HEPGEOMETRYSVC_H
#define ATS_DD4HEPGEOMETRYINTERFACES_IDD4HEPGEOMETRYSVC_H 1

// Gaudi
#include "GaudiKernel/IService.h" 

namespace DD4hep {
    namespace Geometry {
        class DetElement;
    }
}

namespace Add4hep {
    
    static const InterfaceID IID_IDD4hepGeometrySvc("IDD4hepGeometrySvc", 1, 0);
    
    /** @class IDD4hepGeometrySvc
        
        Interface for the service providing the DD4hep geometry.
     
        @author julia.hrdinka@cern.ch
     
     */
    
    class IDD4hepGeometrySvc : virtual public IService {
    
    public:
        
        /** virtual destructor */
        virtual ~IDD4hepGeometrySvc(){}
        
        /** service interface method */
        static const InterfaceID& interfaceID() { return IID_IDD4hepGeometrySvc; }
        
        /** returns the DD4hep world detector element */
        virtual const DD4hep::Geometry::DetElement& worldDetElement() const = 0;
    };
}

#endif //ATS_DD4HEPGEOMETRYINTERFACES_IDD4HEPGEOMETRYSVC_H

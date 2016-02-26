///////////////////////////////////////////////////////////////////
// DD4hepDetectorElement.h, ATS project, DD4hepDetector plugin
///////////////////////////////////////////////////////////////////

#ifndef ATS_DD4HEPDETECTORELEMENT_DD4HEPDETECTORELEMENT_H
#define ATS_DD4HEPDETECTORELEMENT_DD4HEPDETECTORELEMENT_H 1

//TGeoDetectorElement
#include "TGeoDetectorElement/TGeoDetectorElement.h"
//DD4hep
#include "DD4hep/Detector.h"

namespace Ats {
    
    /** @class DD4hepDetectorElement
     
     DetectorElement plugin for DD4hep detector elements. DD4hep is based on TGeo shapes, therefore the DD4hepDetectorElement inherits from TGeoDetectorElement. The full geometrical information is provided by the TGeoDetectorElement. Only the TGeoShape and the DD4hepIdentifier need to be provided.
     
     @author julia.hrdinka@cern.ch
     @TODO what if shape conversion failes? add implementation of more than one surface per module, implementing also for other shapes->Cone,ConeSeg,Tube? what if not used with DD4hep?
     
     */
    
    class DD4hepDetectorElement : public TGeoDetectorElement {
        
    public:
        
        /** Constructor */
        DD4hepDetectorElement(const DD4hep::Geometry::DetElement& detElement, const DD4hep::Geometry::Segmentation& segmentation, std::shared_ptr<const Ats::Transform3D> motherTransform=nullptr);
        /** Desctructor */
        virtual ~DD4hepDetectorElement();
        
    private:
        /** DD4hep detector element */
        DD4hep::Geometry::DetElement            m_detElement;
        /** DD4hep segmentation */
        DD4hep::Geometry::Segmentation          m_segmentation;
        
    };
}


#endif //ATS_DD4HEPDETECTORELEMENT_DD4HEPDETECTORELEMENT_H

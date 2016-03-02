///////////////////////////////////////////////////////////////////
// DetSensComponent.h, ATS project, DD4hepDetector plugin
///////////////////////////////////////////////////////////////////

#ifndef ATS_DD4HEPDETECTORELEMENT_DETSENSCOMPONENT_H
#define ATS_DD4HEPDETECTORELEMENT_DETSENSCOMPONENT_H 1

#include "DD4hepDetectorElement/IDetExtension.h"
#include <memory>

namespace DD4hep {
    namespace Geometry {
        class DetElement;
        class Segmentation;
    }
}

namespace Add4hep {
    
    /** @class DetSensComponent
     
     This class uses the extension mechanism of DD4hep and provides the Segmentation of a Detector Component to the tracking geometry
     
     @author julia.hrdinka@cern.ch
     */
    
    class DetSensComponent : public IDetExtension {
    
    public:
        
        DetSensComponent(const DD4hep::Geometry::Segmentation segmentation) :
        m_segmentation(segmentation)
        {}
        DetSensComponent (const DetSensComponent&, const DD4hep::Geometry::DetElement&)
        {}
        virtual ~DetSensComponent()
        {}
        const DD4hep::Geometry::Segmentation segmentation()
        {
            return (m_segmentation);
        }
        virtual ExtensionType type()
        {
            return ExtensionType::SensComponent;
        }
        
    private:
        
        const DD4hep::Geometry::Segmentation m_segmentation;
    };
}

#endif //ATS_DD4HEPDETECTORELEMENT_DETSENSCOMPONENT_H

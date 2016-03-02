///////////////////////////////////////////////////////////////////
// DetCylinderVolume.h, ATS project, DD4hepDetector plugin
///////////////////////////////////////////////////////////////////

#ifndef ATS_DD4HEPDETECTORELEMENT_DETCYLINDERVOLUME_H
#define ATS_DD4HEPDETECTORELEMENT_DETCYLINDERVOLUME_H

#include "DD4hepDetectorElement/IDetExtension.h"

namespace DD4hep {
    namespace Geometry {
        class DetElement;
    }
}

namespace Add4hep {
    
    /** @class DetDiscVolume
     
     This class uses the extension mechanism of DD4hep to distinguish in the translation to the tracking geometry between a disc and a cylinder volume, since in DD4hep they are both described with the ROOT TGeoConeElement class.
     
     @author julia.hrdinka@cern.ch
     */
    
    class DetCylinderVolume : public IDetExtension {
        
    public:
        
        DetCylinderVolume()
        {}
        DetCylinderVolume(const DetCylinderVolume& volume, const DD4hep::Geometry::DetElement&)
        {}
        virtual ~DetCylinderVolume()
        {}

        virtual ExtensionType type()
        {
            return ExtensionType::CylinderVolume;
        }
        
    };
}

#endif //ATS_DD4HEPDETECTORELEMENT_DETCYLINDERVOLUME_H

///////////////////////////////////////////////////////////////////
// IDetExtension.h, ATS project, DD4hepDetector plugin
///////////////////////////////////////////////////////////////////

#ifndef ATS_DD4HEPDETECTORELEMENT_IDETEXTENSION_H
#define ATS_DD4HEPDETECTORELEMENT_IDETEXTENSION_H 1


namespace Add4hep {
    
    /** @class IDetExtension
     
     Interface class for making extensions to the DD4hep::DetElement class.
     
     @author julia.hrdinka@cern.ch
     */
    
    enum ExtensionType {
        
        None           = 0,
        CylinderVolume = 1,
        DiscVolume     = 2,
        SensComponent  = 3
    };
    
    class IDetExtension {
    
    public:
        
        virtual ~IDetExtension()
        {}
        virtual ExtensionType type() = 0;
        
    protected:
        
        IDetExtension()
        {}
    };
}

#endif //ATS_DD4HEPDETECTORELEMENT_DET_IDETEXTENSION_H

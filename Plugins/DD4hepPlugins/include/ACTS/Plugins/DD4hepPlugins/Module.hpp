///////////////////////////////////////////////////////////////////
// Module.h, ACTS project, Modules plugin
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DD4HEPDETECTORELEMENT_MODULE_H
#define ACTS_DD4HEPDETECTORELEMENT_MODULE_H 1

#include <memory>
//Algebra
#include "ACTS/Utilities/Definitions.hpp"
//DD4hep
#include "DD4hep/Detector.h"


namespace Acts {
    
    /** @class Module
     
    Helper Class needed for the material approxiamtion of a detector module. Stores one DD4hep::Geometry::DetElement, which describes a detector module plus its various placements.
     @TODO find replacement for Gaudi exeption and message stream
     
     @author julia.hrdinka@cern.ch
     */
    
    class Module {
    
    public:
        
        /**constructor*/
        Module(const DD4hep::Geometry::DetElement& mod, std::vector<std::shared_ptr<const Acts::Transform3D>> placements);
        
        /**destructor*/
        ~Module();
        
        /** Assignment operator */
        Module& operator=(const Module& mod);
        
        /**return the module*/
        const DD4hep::Geometry::DetElement& module() const;
        
        /**return the placements*/
        std::vector<std::shared_ptr<const Acts::Transform3D>> placements() const;
        
    private:
        DD4hep::Geometry::DetElement                                 m_mod; //!<detector module
        std::vector<std::shared_ptr<const Acts::Transform3D>>        m_placements; //!<< corresponding placements
    };
}

inline const DD4hep::Geometry::DetElement& Acts::Module::module() const {
    return m_mod;
}

inline std::vector<std::shared_ptr<const Acts::Transform3D>> Acts::Module::placements() const {
    return m_placements;
}

#endif //ACTS_DD4HEPDETECTORELEMENT_MODULE_H

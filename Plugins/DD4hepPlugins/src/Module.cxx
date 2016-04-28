#include "ACTS/Plugins/DD4hepPlugins/Module.h"

Acts::Module::Module(const DD4hep::Geometry::DetElement& mod, std::vector<std::shared_ptr<const Acts::Transform3D>> placements) :
    m_mod(mod),
    m_placements(placements)
    {}

Acts::Module::~Module()
{}

Acts::Module& Acts::Module::operator=(const Acts::Module& mod) {
    
    if (this!=&mod){
        m_mod        = mod.m_mod;
        m_placements = mod.m_placements;
    }
    return *this;
}





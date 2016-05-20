// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Plugins/DD4hepPlugins/Module.hpp"

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





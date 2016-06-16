// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialStepCollection.hpp, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_MATERIAL_MATERIALSTEPCOLLECTION_H
#define ACTS_MATERIAL_MATERIALSTEPCOLLECTION_H

#include "ACTS/Plugins/MaterialPlugins/MaterialStep.hpp"
#include <vector>


namespace Acts {
    
    /**
     @class MaterialStepCollection
     */
    
    class MaterialStepCollection
    {
    public:
        /** constructor */
        MaterialStepCollection();
        
        MaterialStepCollection(std::vector<MaterialStep> steps);
        /** copy constructor */
        MaterialStepCollection(const MaterialStepCollection& msteps);
        /** Pseudo-Constructor clone() */
        MaterialStepCollection* clone() const;
        /** assignment operator */
        MaterialStepCollection& operator=(const MaterialStepCollection& msteps);
        /** destructor */
        ~MaterialStepCollection() = default;
        
        std::vector<MaterialStep> m_steps;
    };
}

#endif //#ifndef ACTS_MATERIAL_MATERIALSTEPCOLLECTION_H
// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialStepCollection.cpp, ACTS project
///////////////////////////////////////////////////////////////////
#include "ACTS/Plugins/MaterialPlugins/MaterialStepCollection.hpp"

Acts::MaterialStepCollection::MaterialStepCollection() :
m_steps()
{}

Acts::MaterialStepCollection::MaterialStepCollection(std::vector<MaterialStep> steps) :
m_steps(steps)
{}


Acts::MaterialStepCollection::MaterialStepCollection(const MaterialStepCollection& msteps) :
m_steps(msteps.m_steps)
{}

Acts::MaterialStepCollection* Acts::MaterialStepCollection::clone() const
{
    return (new MaterialStepCollection(*this));
}

Acts::MaterialStepCollection& Acts::MaterialStepCollection::operator=(const MaterialStepCollection& msteps)
{
    if(this != &msteps)
    {
        m_steps     = msteps.m_steps;
    }
    return *this;
}

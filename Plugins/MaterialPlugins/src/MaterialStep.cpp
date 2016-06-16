// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialStep.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Plugins/MaterialPlugins/MaterialStep.hpp"

Acts::MaterialStep::MaterialStep() :
X0(0.),
L0(0.),
A(0.),
Z(0.),
rho(0.),
position(Position()),
steplength(0.)
{}

Acts::MaterialStep::MaterialStep(double X0, double L0, double A, double Z, double rho, Position pos, double steplength) :
X0(X0),
L0(L0),
A(A),
Z(Z),
rho(rho),
position(pos),
steplength(steplength)
{}

Acts::MaterialStep::MaterialStep(const MaterialStep& mstep) :
X0(mstep.X0),
L0(mstep.L0),
A(mstep.A),
Z(mstep.Z),
rho(mstep.rho),
position(mstep.position),
steplength(mstep.steplength)
{}

Acts::MaterialStep* Acts::MaterialStep::clone() const
{
    return (new MaterialStep(*this));
}

Acts::MaterialStep& Acts::MaterialStep::operator=(const MaterialStep& mstep)
{
    if(this != &mstep)
    {
        X0              = mstep.X0;
        L0              = mstep.L0;
        A               = mstep.A;
        Z               = mstep.Z;
        rho             = mstep.rho;
        position        = mstep.position;
        steplength      = mstep.steplength;
    }
    return *this;
}

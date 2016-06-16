// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialStep.hpp, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_MATERIAL_MATERIALSTEP_H
#define ACTS_MATERIAL_MATERIALSTEP_H

#include <memory>


namespace Acts {
    
    /**
     @class MaterialStep
     */
    
    class MaterialStep
    {
    public:
        
        struct Position {
            double x;
            double y;
            double z;
            
            Position() :
            x(0.),
            y(0.),
            z(0.)
            {}
            
            Position(double x, double y, double z) :
            x(x),
            y(y),
            z(z)
            {}
        };
        /** constructor */
        MaterialStep();
        
        MaterialStep(double X0, double L0, double A, double Z, double rho, Position pos, double steplength);
        /** copy constructor */
        MaterialStep(const MaterialStep& mstep);
        /** Pseudo-Constructor clone() */
        MaterialStep* clone() const;
        /** assignment operator */
        MaterialStep& operator=(const MaterialStep& mstep);
        /** destructor */
        ~MaterialStep() = default;
        
        
        double                                          X0; //!<< material properties
        double                                          L0;
        double                                          A;
        double                                          Z;
        double                                          rho;
        Position                                        position;
        double                                          steplength;
        
    };
}

#endif //ACTS_MATERIAL_MATERIALSTEP_H
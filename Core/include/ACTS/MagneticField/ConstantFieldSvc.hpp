// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ConstantFieldSvc.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_MAGNETICFIELDINTERFACES_CONSTANTFIELDSVC_H
#define ACTS_MAGNETICFIELDINTERFACES_CONSTANTFIELDSVC_H 1

#include "ACTS/MagneticField/IMagneticFieldSvc.hpp"
#include <array>
#include <string>

namespace FWE {
    
    /// @class ConstantFieldSvc
    ///
    /// @brief returns a given constant field value at every point
    ///
    ///  Implementation of the IMagneticFieldSvc interface with a constant magnetic field.
    ///  The Value can be set at configuration in kilo Tesla. This service returns the
    ///  constant magnetic field at any point.
    ///
    
    class ConstantFieldSvc : public Acts::IMagneticFieldSvc {
        
        ///////////////////////////////////////////////////////////////////
        // Public methods:
        ///////////////////////////////////////////////////////////////////
    public:
        /// @class Config - nested configuraiton class
        class Config {
        public:
            std::array<double,3> field;
            std::string          name;
            
            Config():
             field({{0.,0.,20.}}),
             name("Anonymous")
            {}
        };
        
        /// Constructor
        ConstantFieldSvc(const Config& cfg):
         m_cfg(cfg)
        {}
        
        /// Destructor
        ~ConstantFieldSvc(){}
        
        /// get B field value at given position
        /// @param xyz xyz[3] is the position in mm,
        /// @param bxyz bxyz[3] is the magnetic field value in kT
        /// @param deriv if deriv[9] is given, field derivatives are returned in kT/mm
        void getField(const double *xyz, double *bxyz, double *deriv = nullptr) const final;
        
        /// get B field value on the z-r plane at given position */
        /// @note works only inside the solenoid; otherwise calls getField() above */
        /// @param xyz xyz[3] is the position in mm,
        /// @param bxyz bxyz[3] is the magnetic field value in kT
        /// @param deriv if deriv[9] is given, field derivatives are returned in kT/mm
        void getFieldZR(const double *xyz, double *bxyz, double *deriv = nullptr) const final;
        
    private:
        Config m_cfg;
        
    };
    
    inline void ConstantFieldSvc::getField(const double*, double* bxyz, double *) const
    {
        bxyz[0] = m_cfg.field[0];
        bxyz[1] = m_cfg.field[1];
        bxyz[2] = m_cfg.field[2];
        return;
    }

    inline void ConstantFieldSvc::getFieldZR(const double*, double* bxyz, double *) const
    {
        bxyz[0] = m_cfg.field[0];
        bxyz[1] = m_cfg.field[1];
        bxyz[2] = m_cfg.field[2];
        return;
    }
    
    
}


#endif //> ! ACTS_MAGNETICFIELDINTERFACES_CONSTANTFIELDSVC_H


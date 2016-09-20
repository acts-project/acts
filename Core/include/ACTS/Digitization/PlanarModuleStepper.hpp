// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_DIGITIZATION_PLANARMODULESTEPPER_H
#define ACTS_DIGITIZATION_PLANARMODULESTEPPER_H

#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Digitization/DigitizationCell.hpp"
#include <vector>
#include <memory>

namespace Acts {

    class DigitizationModule;
    
    /// @class PlanarModuleStepper
    ///
    /// Module for fast, geometric digitization
    /// this is a planar module stepper that calculates the step length
    /// in given segmentations and retrunes digitisation steps

    class PlanarModuleStepper {

      public:
        ///  @struct Config
        ///  Configuration for the planar module stepper
        struct Config
        {
          // standard constructor
          Config(){}
        };
      
        /// Constructor
        /// @param pmsConfig is the configuration lo
        PlanarModuleStepper(const Config&  pmsConfig,
                     std::unique_ptr<Logger> logger
                     = getDefaultLogger("LayerCreator", Logging::INFO));
      
        /// Destructor
        ~PlanarModuleStepper() = default;

        /// calculate the steps caused by this track - full simulation interface 
        std::vector<DigitizationStep> cellSteps(const DigitizationModule& dmodule,
                                                const Vector3D& startPosition,
                                                const Vector3D& endPosition) const;
        
        /// calculate the steps caused by this track - fast simulation interface */
        std::vector<DigitizationStep> cellSteps(const DigitizationModule& dmodule,
                                                const Vector2D& intersection,
                                                const Vector3D& direction) const;
                                                
        /// set logging instance
        void
        setLogger(std::unique_ptr<Logger> logger)
        {
          m_logger = std::move(logger);
        }
     
      private:
        /// Private access method to the logging instance
        const Logger&
        logger() const
        {
          return *m_logger;
        }
     
        /// logging instance
        std::unique_ptr<Logger> m_logger;

    };

} // end of namespace

#endif // ACTS_DIGITIZATION_PLANARMODULESTEPPER_H

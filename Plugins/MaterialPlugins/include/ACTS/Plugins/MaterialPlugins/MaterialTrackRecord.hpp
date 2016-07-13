// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialTrackRecord.hpp, ACTS project MaterialPlugins
///////////////////////////////////////////////////////////////////

#ifndef ACTS_MATERIAL_MATERIALTRACKRECORD_H
#define ACTS_MATERIAL_MATERIALTRACKRECORD_H

#include "ACTS/Plugins/MaterialPlugins/MaterialStep.hpp"
#include <vector>


namespace Acts {
    
    
    /// @class MaterialTrackRecord
    ///
    /// @brief holds the material steps along a track, the direction and the starting point of the track
    ///
    /// The class MaterialTrackRecord holds a collection of material steps (containg the material in one step and the
    /// position) along a track through the detector, the three dimensional global starting point of the track and the
    /// direction expressed in pseudorapidity eta and azimuthal angle phi.
    ///
    /// For this class a ROOT dictionary is created in order to store it in a ROOT tree.
    
    class MaterialTrackRecord {
    
    public:
        /// default constructor
        MaterialTrackRecord() = default;
        /// Constructor handing over all the relevant parameters for the material track
        /// @param startPos three dimensional global start position of the track
        /// @param eta pseudorapidity indicating the first coordinate of the direction of the material track
        /// @param phi azimuthal angle indicating the second coordinate of the direction of the material track
        /// @param materialSteps the collection material steps along the track
        MaterialTrackRecord(const MaterialStep::Position& startPos, double eta, double phi, std::vector< MaterialStep> materialSteps);
        /// Copy constructor
        MaterialTrackRecord(const MaterialTrackRecord& mtrecord);
        /// Implicit contructor
        /// - uses the copy constructor
        const MaterialTrackRecord* clone() const;
        /// Default destructor
        ~MaterialTrackRecord() = default;
        /// Assignment operator
        MaterialTrackRecord& operator=(const MaterialTrackRecord& mtrecord);
        /// Return method for the pseudorapidity eta
        double eta() const;
        /// Return method for the azimuthal angle phi
        double phi() const;
        /// Return method for the start position of the track
        const MaterialStep::Position position() const;
        /// Return method for the material step
        std::vector<MaterialStep> materialSteps() const;
        
    private:
        /// start position of the material track
        MaterialStep::Position                m_startPosition;
        /// pseudorapidity indicating the first coordinate of the direction of the material track
        double                                      m_eta;
        /// azimuthal angle phi indicating the second coordinate of the direction of the material track
        double                                      m_phi;
        /// the collection material steps along the track
        std::vector<MaterialStep>                   m_materialSteps;
    };
}

inline double Acts::MaterialTrackRecord::eta() const
{
    return m_eta;
}

inline double Acts::MaterialTrackRecord::phi() const
{
    return m_phi;
}

inline const Acts::MaterialStep::Position Acts::MaterialTrackRecord::position() const
{
    return m_startPosition;
}

inline std::vector<Acts::MaterialStep> Acts::MaterialTrackRecord::materialSteps() const
{
    return m_materialSteps;
}

#endif //ACTS_MATERIAL_MATERIALTRACKRECORD_H
// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DD4hepCylinderGeometryBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H
#define ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H 1

// Geometry module
#include "ACTS/Tools/ITrackingGeometryBuilder.hpp"
#include "ACTS/Tools/ITrackingVolumeBuilder.hpp"
#include "ACTS/Tools/ITrackingVolumeHelper.hpp"
#include "ACTS/Layers/Layer.hpp"
// DD4hep
#include "DD4hep/Detector.h"

namespace Acts {
    class TrackingGeometry;
}
namespace Acts {

    
    /** @ class DD4hepCylinderGeometryBuilder
     
    This class receives the DD4hep geometry from the given implementation of the IDD4hepGeomtrySvc, walks through the subvolumes and initiates their translation into the tracking geometry.
        It returns the world tracking geometry element.
     @TODO find replacement for Gaudi exeption and message stream
     
     */
    
    
    class DD4hepCylinderGeometryBuilder : virtual public Acts::ITrackingGeometryBuilder {
        
    public:
        /** @struct Config
         Configuration for the DD4hepCylinderGeometryBuilder */
        struct Config {
        
            std::shared_ptr<ITrackingVolumeBuilder>     volumeBuilder; //!< building the contained sub detectors
            std::shared_ptr<ITrackingVolumeHelper>      volumeHelper; //!< helper tool needed for volume building
            DD4hep::Geometry::DetElement                detWorld;     //!< world detector element of the DD4hep geometry
            
            Config() :
                volumeBuilder(nullptr),
                volumeHelper(nullptr),
                detWorld()
            {}
        };
        
        /** Constructor */
        DD4hepCylinderGeometryBuilder(const Config dgbConfig);
        
        /** Destructor */
        virtual ~DD4hepCylinderGeometryBuilder();
        
        /** setting the builders and helpers with the configuration object*/
        void setConfiguration(const Config dgbConfig);
        
        /** get the configuration object */
        Config getConfiguration() const;
        
        /** TrackingGeometry Interface method - translates the DD4hep geometry into the tracking geometry*/
        std::unique_ptr<TrackingGeometry> trackingGeometry() const override;
        
    private:
        /** configuration object */
        Config                                            m_config;
        
    };
    
    inline DD4hepCylinderGeometryBuilder::Config DD4hepCylinderGeometryBuilder::getConfiguration() const
    { return m_config; }
} //end of namespace

#endif //#ifndef ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H

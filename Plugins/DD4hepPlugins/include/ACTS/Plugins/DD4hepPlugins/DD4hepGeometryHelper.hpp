// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DD4hepGeometryHelper.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPGEOMETRYHELPER_H
#define ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPGEOMETRYHELPER_H

// Core module
#include "ACTS/Utilities/Definitions.hpp"
// Geometry Module
#include "ACTS/Volumes/VolumeBounds.hpp"
// Geometry module
#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Tools/ITrackingVolumeBuilder.hpp"
// DD4hep
#include "DD4hep/Detector.hpp"


namespace Acts {
    
    /** @ class DD4hepGeometryHelper
     
     Provides helper function to translate the DD4hep geometry into the ACTS Tracking Geometry.
     @TODO find replacement for Gaudi exeption and message stream
     
     */

   class DD4hepGeometryHelper  {
       
   public:
       /** constructor */
       DD4hepGeometryHelper() = default;
       
       /** destructor */
       ~DD4hepGeometryHelper() = default;
       
       /**helper method to extract the transformation matrix from a DD4hep DetElement*/
       static std::shared_ptr<Acts::Transform3D> extractTransform(DD4hep::Geometry::DetElement& detElement);
       /**helper method to extract the volume boundaries of a cylindrical volume*/
       static std::shared_ptr<const Acts::VolumeBounds> extractVolumeBounds(DD4hep::Geometry::DetElement& detElement);
       /** Creates a triple of volumes a possible barrel-endcap configuration and of all the three possible Layer types of the given volume detector element*/
       /** constructs all subvolumes contained by this volume (motherDetELement) with its layers and modules, if present */
       static void createSubVolumes(DD4hep::Geometry::DetElement& motherDetElement, LayerTriple& layerTriple, VolumeTriple& volumeTriple);
       
   private:
       /**creates the cylindrical shaped layers*/
       static void createCylinderLayers(DD4hep::Geometry::DetElement& motherDetElement, Acts::LayerVector& centralLayers, std::shared_ptr<Acts::Transform3D> motherTransform = nullptr);
       /**creates disc shaped layers*/
       static void createDiscLayers(DD4hep::Geometry::DetElement& motherDetElement, Acts::LayerVector& layers, std::shared_ptr< Acts::Transform3D> motherTransform = nullptr);
       /**creates a binned array of Acts::Surfaces out of vector of DD4hep detector modules*/
       static std::unique_ptr<Acts::SurfaceArray> createSurfaceArray(std::vector<DD4hep::Geometry::DetElement>& modules, Acts::BinningValue lValue, std::shared_ptr<const Acts::Transform3D> motherTransform = nullptr);
       /**creating a surface array binned in phi and a longitudinal direction which can either be z or r*/
       static std::unique_ptr<Acts::SurfaceArray> binnedSurfaceArray2DPhiL(const std::vector<const Acts::Surface*> surfaces, Acts::BinningValue lValue);
       /**helper method to get the bin values for a binned array out of overlapping modules*/
       static std::vector<float> createBinValues(std::vector<std::pair<float,float>> old);
        /**helper method to sort pairs of doubles*/
       static bool sortFloatPairs(std::pair<float,float> ap, std::pair<float,float> bp);
       
    };
}

#endif //ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPGEOMETRYHELPER_H

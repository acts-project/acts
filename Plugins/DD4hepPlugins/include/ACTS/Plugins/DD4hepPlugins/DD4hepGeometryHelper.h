///////////////////////////////////////////////////////////////////
// DD4hepLayerHelper.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPLAYERHELPER_H
#define ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPLAYERHELPER_H 1

// Core module
#include "ACTS/Utilities/Definitions.h"
// Geometry module
#include "ACTS/Detector/TrackingVolume.h"
#include "ACTS/Utilities/BinnedArray.h"
#include "ACTS/Tools/ITrackingVolumeBuilder.h"
#include "ACTS/Plugins/DD4hepPlugins/Module.h"
// DD4hep
#include "DD4hep/Detector.h"

namespace Acts {
    class Surface;
    class Volume;
}

namespace Acts {
    
    /** @ class DD4hepGeometryHelper
     
     The DD4hepGeometryHelper allows to create the subvolumes in the tracking geometry of a barrel or barrel-endcap hierarchy given from DD4hep. If present it creates a volume triple <nEndcap,Barrel,pEndcap> and if present its corresponing layers with possibly containing Modules and/or support structure also as a triple. The triples contain pointers to the particular members, which are set to zero if a certain member is not present.

     @TODO find replacement for Gaudi exeption and message stream
     
     @author julia.hrdinka@cern.ch
    */
    
    typedef std::vector<const Acts::Surface*>                        SurfaceVector;
    typedef Acts::BinnedArray<const Acts::Surface*>                   SurfaceArray;
    
   class DD4hepGeometryHelper  {
       
   public:
       /** constructor */
       DD4hepGeometryHelper();
       
       /** destructor */
       ~DD4hepGeometryHelper();
       
       /**helper method to extract the transformation matrix from a DD4hep DetElement*/
       static std::shared_ptr<Acts::Transform3D> extractTransform(DD4hep::Geometry::DetElement& detElement);
       /**helper method to extract the volume boundaries of a cylindrical volume*/
       static std::shared_ptr<const Acts::VolumeBounds> extractVolumeBounds(DD4hep::Geometry::DetElement& detElement);
       /** Creates a triple of volumes a possible barrel-endcap configuration and of all the three possible Layer types of the given volume detector element*/
       /** constructs all subvolumes contained by this volume (motherDetELement) with its layers and modules, if present */
       void createSubVolumes(DD4hep::Geometry::DetElement& motherDetElement, LayerTriple& layerTriple, VolumeTriple& volumeTriple);
       
   private:
       /**creates the cylindrical shaped layers*/
       void createCylinderLayers(DD4hep::Geometry::DetElement& motherDetElement, Acts::LayerVector& centralLayers, std::shared_ptr<Acts::Transform3D> motherTransform = nullptr);
       /**creates disc shaped layers*/
       void createDiscLayers(DD4hep::Geometry::DetElement& motherDetElement, Acts::LayerVector& layers, std::shared_ptr< Acts::Transform3D> motherTransform = nullptr);
       /**creates a binned array of Acts::Surfaces out of vector of DD4hep detector modules*/
       createSurfaceArray(std::vector<Module>& modules, Acts::BinningValue lValue, std::shared_ptr<const Acts::Transform3D> motherTransform = nullptr) const
       /**creating a surface array binned in phi and a longitudinal direction which can either be z or r*/
       binnedSurfaceArray2DPhiL(const std::vector<const Acts::Surface*>& surfaces, Acts::BinningValue lValue) const
       /**helper method to get the bin values for a binned array out of overlapping modules*/
       std::vector<float> createBinValues(std::vector<std::pair<float,float>> old) const;
        /**helper method to sort pairs of doubles*/
       static bool sortFloatPairs(std::pair<float,float> ap, std::pair<float,float> bp);
       
    };
    
} //end of namespace Acts


#endif //ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPGEOMETRYHELPER_H

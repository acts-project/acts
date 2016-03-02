///////////////////////////////////////////////////////////////////
// DD4hepLayerHelper.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_DD4HEPGEOMETRYTOOLS_DD4HEPLAYERHELPER_H
#define ATS_DD4HEPGEOMETRYTOOLS_DD4HEPLAYERHELPER_H 1

//Interface
#include "DD4hepGeometryInterfaces/IDD4hepLayerHelper.h"
// Core module
#include "CoreInterfaces/AlgToolBase.h"
#include "Algebra/AlgebraDefinitions.h"
// Geometry module
#include "Detector/TrackingVolume.h"
#include "GeometryUtils/BinnedArray.h"
// Gaudi & Athena
#include "GaudiKernel/ToolHandle.h"
// DD4hep
#include "DD4hep/Detector.h"

namespace Ats {
    class Surface;
}

namespace Add4hep {
    
    /** @ class DD4hepLayerHelper
        
     The DD4hepLayerHelper translates Layers with its containing modules (if given) from the DD4hep geometry into Layers and Surfaces from TrackingGeometry.
     
     @author julia.hrdinka@cern.ch
    */
    
    class DD4hepGeometryHelper;
    
    typedef std::vector<const Ats::Surface*>                        SurfaceVector;
    typedef Ats::BinnedArray<const Ats::Surface*>                   SurfaceArray;
    
   class DD4hepLayerHelper : public Ats::AlgToolBase, virtual public IDD4hepLayerHelper {
       
   public:
       /** constructor */
       DD4hepLayerHelper(const std::string& t, const std::string& n, const IInterface* p);
       
       /** destructor */
       ~DD4hepLayerHelper();
       
       /** AlgTool initilaize method */
       StatusCode initialize() override;
       
       /** AlgTool finalize method*/
       StatusCode finalize() override;
       
       /** Creates a triple of all the three possible Layer types of the given volume detector element*/
       virtual const Ats::LayerTriple* createLayerTriple(DD4hep::Geometry::DetElement& motherDetElement) override;
       
   private:
       /**constructs all layers contained by this volume*/
       StatusCode constructLayers(DD4hep::Geometry::DetElement& detElement);
       /**creates the cylindrical shaped layers*/
       StatusCode createCylinderLayers(DD4hep::Geometry::DetElement& motherDetElement, std::shared_ptr<const Ats::Transform3D> motherTransform) const;
       /**creates disc shaped layers*/
       StatusCode createDiscLayers(DD4hep::Geometry::DetElement& motherDetElement, std::shared_ptr<const Ats::Transform3D> motherTransform) const;
       /**method to create a binned array out of  a vector of surfaces for cylinder layers (surfaces binned in z and phi)*/
       SurfaceArray* createCylinderBinnedSurfaceArray(DD4hep::Geometry::DetElement& motherDetElement, std::shared_ptr<const Ats::Transform3D> motherTransform, double Zmin, double Zmax) const;
       /**method to create a binned array out of a vector of surfaces for disc layers (surfaces binned in r and phi)*/
       SurfaceArray* createDiscBinnedSurfaceArray(DD4hep::Geometry::DetElement& motherDetElement, std::shared_ptr<const Ats::Transform3D> motherTransform) const;
       /**method which creates a vector of contained surfaces by a given DD4hep mother DetElement - used internally to extract the surfaces of a layer*/
       StatusCode createSurfaceVector(DD4hep::Geometry::DetElement& motherDetElement, std::shared_ptr<const Ats::Transform3D> motherTransform, SurfaceVector& surfaces) const;
       /**helper method to get the bin values for a binned array out of overlapping modules*/
       std::vector<float> orderRValues(std::vector<std::pair<float,float>> old) const;
        /**helper method to sort pairs of doubles*/
       static bool sortFloatPairs(std::pair<float,float> ap, std::pair<float,float> bp);
       
       /**possible layers of the negative end cap*/
       Ats::LayerVector*                                m_negativeLayers;
       /**possible layers of the central barrel*/
       Ats::LayerVector*                                m_centralLayers;
       /**possible layers of the positive end cap*/
       Ats::LayerVector*                                m_positiveLayers;
       
       
    };
    
} //end of namespace Add4hep



#endif //ATS_DD4HEPGEOMETRYTOOLS_DD4HEPLAYERHELPER_H

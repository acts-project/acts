///////////////////////////////////////////////////////////////////
// CylinderVolumeHelper.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYTOOLS_CYLINDERVOLUMEHELPER_H
#define ACTS_GEOMETRYTOOLS_CYLINDERVOLUMEHELPER_H

#ifndef TRKDETDESCR_TAKESMALLERBIGGER
#define TRKDETDESCR_TAKESMALLERBIGGER
#define takeSmaller(current,test) current = current < test ? current : test
#define takeBigger(current,test)  current = current > test ? current : test
#define takeSmallerBigger(cSmallest, cBiggest, test) takeSmaller(cSmallest, test); takeBigger(cBiggest, test)
#endif

// Core module
#include "CoreInterfaces/AlgToolBase.h"
// Geometry module
#include "GeometryInterfaces/ITrackingVolumeHelper.h"
#include "Volumes/BoundarySurfaceFace.h"
#include "GeometryInterfaces/ILayerArrayCreator.h"
#include "GeometryInterfaces/ITrackingVolumeArrayCreator.h"
// Gaudi & Athena
#include "GaudiKernel/ToolHandle.h"
// STL
#include <string>
#include <vector>

namespace Acts {
    
    class Layer;
    class TrackingVolume;
    class VolumeBounds;
    class CylinderVolumeBounds;
    class Material;
    
    /** @class CylinderVolumeHelper
     
     The concrete implementation for cylindrical TrackingVolume
     objects of the ITrackingVolumeCreator interface
     
     @TODO Julia: make private tools private again after Gaudi update (bug in Gaudi), marked with //b
     
     @author Andreas.Salzburger@cern.ch
     */
    
    class CylinderVolumeHelper : public AlgToolBase, virtual public ITrackingVolumeHelper {
        
    public:
        /** Constructor */
        CylinderVolumeHelper(const std::string&,const std::string&,const IInterface*);
        
        /** Destructor */
        ~CylinderVolumeHelper();
        
        /** AlgTool initialize method */
        virtual StatusCode initialize() override;
        
        /** AlgTool finalize method */
        virtual StatusCode finalize() override;
        
        /** @copydoc ITrackingVolumeCreator::createTrackingVolume(const LayerVector&, const Material& matprop, VolumeBounds*, Transform3D*,bool, const std::string&) const;  */
        TrackingVolumePtr createTrackingVolume(const LayerVector& layers,
                                               const Material& matprop,
                                               VolumeBoundsPtr volBounds,
                                               std::shared_ptr<Transform3D> transform = nullptr,
                                               const std::string& volumeName = "UndefinedVolume",
                                               BinningType btype = arbitrary) const;
        
        /** @copydoc ITrackingVolumeCreator::createTrackingVolume(const std::vector<const Layer*>& , const Material&, ,double,double,double,double,bool,const std::string&) const;
         */
        TrackingVolumePtr createTrackingVolume(const LayerVector& layers,
                                               const Material& matprop,
                                               double loc1Min, double loc1Max,
                                               double loc2Min, double loc2Max,
                                               const std::string& volumeName = "UndefinedVolume",
                                               BinningType btype = arbitrary) const;
        
        /** @copydoc ITrackingVolumeCreator::createGapTrackingVolume(const Material&, double,double,double,double,int,bool,const std::string&) const; */
        TrackingVolumePtr createGapTrackingVolume(const Material& matprop,
                                                  double rMin, double rMax,
                                                  double zMin, double zMax,
                                                  unsigned int materialLayers,
                                                  bool cylinder = true,
                                                  const std::string& volumeName = "UndefinedVolume") const;
        
        /** @copydoc ITrackingVolumeCreator::createGaoTrackingVolume(Material&,,std::vector<double>&,int,bool,const std::string&) const;  */
        TrackingVolumePtr createGapTrackingVolume(const Material& matprop,
                                                  double rMin, double rMax,
                                                  double zMin, double zMax,
                                                  const std::vector<double>& layerPositions,
                                                  bool cylinder = true,
                                                  const std::string& volumeName = "UndefinedVolume",
                                                  BinningType btype = arbitrary) const;
        
        /** Create a container volumes from sub volumes, input volumes are ordered in R or Z by convention */
        TrackingVolumePtr createContainerTrackingVolume(const TrackingVolumeVector& volumes) const;
        
        
    private:
        /** Private method - it estimates the CylinderBounds and Translation of layers,
         if given, these are checked against the layer positions/dimensions. */
        StatusCode estimateAndCheckDimension(const LayerVector& layers,
                                             const Acts::CylinderVolumeBounds*& cylBounds,
                                             std::shared_ptr<Transform3D>& transform,
                                             double& rMinClean, double& rMaxClean,
                                             double& zMinClean, double& zMaxClean,
                                             BinningValue& bValue,
                                             BinningType bType = arbitrary) const;
        
        /** Private method - interglue all volumes contained by a TrackingVolume
         and set the outside glue volumes in the descriptor */
        StatusCode interGlueTrackingVolume(TrackingVolumePtr tVolume,
                                           bool rBinned,
                                           double rMin, double rMax,
                                           double zMin, double zMax) const;
        
        /** Private method - glue volume to the other -- use trackingVolume helper */
        void glueTrackingVolumes(TrackingVolumePtr volumeOne,
                                 BoundarySurfaceFace faceOne,
                                 TrackingVolumePtr volumeTwo,
                                 BoundarySurfaceFace faceTwod,
                                 double rMin, double rMax,
                                 double zMin, double zMax) const;
        
        
        /** Private method - helper method not to duplicate code */
        void addFaceVolumes(TrackingVolumePtr tVolume,
                            Acts::BoundarySurfaceFace bsf,
                            TrackingVolumeVector& vols) const;
        
        
        /** Private method - helper method to save some code */
        LayerPtr createCylinderLayer(double z,
                                     double r,
                                     double halflength,
                                     double thickness,
                                     int binsPhi,
                                     int binsZ) const;
        
        /** Private method - helper method to save some code */
        LayerPtr createDiscLayer(double z,
                                 double rMin, double rMax,
                                 double thickness,
                                 int binsPhi,
                                 int binsR) const;
        
        ToolHandle<ILayerArrayCreator>          m_layerArrayCreator;           //!< A Tool for coherent LayerArray creation
        ToolHandle<ITrackingVolumeArrayCreator> m_trackingVolumeArrayCreator;  //!< Helper Tool to create TrackingVolume Arrays
        
        double                                  m_passiveLayerThickness;  //!< thickness of passive layers
        int                                     m_passiveLayerPhiBins;    //!< bins in phi for the passive layer
        int                                     m_passiveLayerRzBins;     //!< bins in r/z for the passive layer
        
    };
}

#endif






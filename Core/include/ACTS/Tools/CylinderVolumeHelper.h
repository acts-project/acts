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

// Geometry module
#include "ACTS/Tools/ITrackingVolumeHelper.h"
#include "ACTS/Volumes/BoundarySurfaceFace.h"
#include "ACTS/Tools/ILayerArrayCreator.h"
#include "ACTS/Tools/ITrackingVolumeArrayCreator.h"
// STL
#include <string>
#include <vector>
#include <memory>

namespace Acts
{    
  class Layer;
  class TrackingVolume;
  class VolumeBounds;
  class CylinderVolumeBounds;
  class Material;
    
  /** @class CylinderVolumeHelper
     
      The concrete implementation for cylindrical TrackingVolume
      objects of the ITrackingVolumeCreator interface
          
      @author Andreas.Salzburger@cern.ch
  */
    
  class CylinderVolumeHelper : public ITrackingVolumeHelper
  {    
  public:
    /** @struct Config 
      Configuration struct for this CylinderVolumeHelper */
    struct Config {  
      std::shared_ptr<ILayerArrayCreator>          layerArrayCreator;           //!< A Tool for coherent LayerArray creation
      std::shared_ptr<ITrackingVolumeArrayCreator> trackingVolumeArrayCreator;  //!< Helper Tool to create TrackingVolume Arrays
      double                                       passiveLayerThickness;       //!< thickness of passive layers
      int                                          passiveLayerPhiBins;         //!< bins in phi for the passive layer
      int                                          passiveLayerRzBins;          //!< bins in r/z for the passive layer    

      Config() :
        layerArrayCreator(nullptr),       
        trackingVolumeArrayCreator(nullptr),
        passiveLayerThickness(1),
        passiveLayerPhiBins(1),
        passiveLayerRzBins(100)
        {}
    };
      
    /** Constructor */
    CylinderVolumeHelper(const Config& cvhConfig);
        
    /** Destructor */
    virtual ~CylinderVolumeHelper() = default;
        
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

    /** Set configuration method */
    void setConfiguration(const Config& cvbConfig);
   
    /** Get configuration method */
    Config getConfiguration() const;   

  protected:
    /** Configuration object */
    Config m_config;  
        
  private:
    /** Private method - it estimates the CylinderBounds and Translation of layers,
	if given, these are checked against the layer positions/dimensions. */
    bool estimateAndCheckDimension(const LayerVector& layers,
				   const Acts::CylinderVolumeBounds*& cylBounds,
				   std::shared_ptr<Transform3D>& transform,
				   double& rMinClean, double& rMaxClean,
				   double& zMinClean, double& zMaxClean,
				   BinningValue& bValue,
				   BinningType bType = arbitrary) const;
        
    /** Private method - interglue all volumes contained by a TrackingVolume
	and set the outside glue volumes in the descriptor */
    bool interGlueTrackingVolume(TrackingVolumePtr tVolume,
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
        
      
  };
  
  /** Return the configuration object */    
  inline CylinderVolumeHelper::Config CylinderVolumeHelper::getConfiguration() const { return m_config; }
}

#endif


///////////////////////////////////////////////////////////////////
// CylinderVolumeBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYTOOLS_CYLINDERVOLUMEBUILDER_H
#define ACTS_GEOMETRYTOOLS_CYLINDERVOLUMEBUILDER_H 1

// Geometry module
#include "ACTS/Tools/ITrackingVolumeBuilder.h"
#include "ACTS/Material/Material.h"
#include "ACTS/Tools/ITrackingVolumeHelper.h"
#include "ACTS/Tools/ILayerBuilder.h"
#include "ACTS/Utilities/BinningType.h"

#ifndef ATAS_GEOMETRYTOOLS_TAKESMALLERBIGGER
#define ATAS_GEOMETRYTOOLS_TAKESMALLERBIGGER
#define takeSmaller(current,test) current = current < test ? current : test
#define takeBigger(current,test)  current = current > test ? current : test
#define takeSmallerBigger(cSmallest, cBiggest, test) takeSmaller(cSmallest, test); takeBigger(cBiggest, test)
#endif

namespace Acts {

  class TrackingVolume;
  class VolumeBounds;
  
  /** @ LayerSetup struct to understand the layer setup */
  struct LayerSetup {      

    bool          present;                        //!< layers are present
    BinningValue  binningValue;                   //!< in what way they are binned
      
    std::pair<double,double> rBoundaries;         //!< raidal boundaries
    std::pair<double,double> zBoundaries;         //!< zBoundaries 
      
    //std::vector<double>      ringBoundaries;      //!< ring boundaries if present //!< @TODO insert ring layout      
    LayerSetup() :
      present(false),
      binningValue(binR),
      rBoundaries(std::pair<double,double>(10e10,-10e10)),
      zBoundaries(std::pair<double,double>(10e10,-10e10))
      {}    
        
    /** Conversion operator to bool */
    operator bool() const { return present; }          
      
  };
    
  /** @class CylinderVolumeBuilder
  
      A simple cylindrical volume builder to be used for building a concentrical cylindrical volume
      - a) configured volume
      - b) wrapping around a cylindrical/disk layer setup
    
      All are optionally wrapped around a given volume which has to by a cylinder volume
      and which has to be center at z == 0
   
      To receive the tracking volume it is possible to also hand over a triple of layers, which is a C++ tuple of three pointers to layer vectors (defined in the ITrackingVolumeBuilder). This functionality is needed for a possible translation of an geometry existing in another format. The first entry represents the layers of the negative endcap, the second the layers of the barrel and the third the layers of the positive endcap. If the one of these pointers is a nullptr no layers will be created for this volume
      Another functionality needed to translate an already existing geometry is to hand over a volume triple, which is a triple of shared pointers of volumes (defined in the ITrackingVolumeBuilder). The first entry contains the negative endcap volume, the second the barrel volume and the third one the positive endcap volume. This volumes are then used to get the internal boundaries of the current hierarchy.
    
      @TODO Julia: make private tools private again after Gaudi update (bug in Gaudi), marked with //b
   
      @author Andreas.Salzburger@cern.ch   
  */
  
  class CylinderVolumeBuilder : public ITrackingVolumeBuilder
  {

  public:
    /** Constructor */
    CylinderVolumeBuilder(const std::string& name);
      
    /** Destructor */
    virtual ~CylinderVolumeBuilder();
      
    /** CylinderVolumeBuilder interface method - returns the volumes of Volumes */
    TrackingVolumePtr trackingVolume(TrackingVolumePtr insideVolume = nullptr,
				     VolumeBoundsPtr outsideBounds  = nullptr,
				     const LayerTriple* layerTriple  = nullptr,
				     const VolumeTriple* volumeTriple = nullptr) const override;

    void setVolumeHelper(std::shared_ptr<ITrackingVolumeHelper> volumeHelper)
    {
      m_trackingVolumeHelper = std::move(volumeHelper);
    }
    
    void setVolumeName(std::string name)
    {
      m_volumeName = std::move(name);
    }

    void setVolumeDimensions(std::vector<double> dim)
    {
      m_volumeDimension = std::move(dim);
    }

    void setVolumeMaterial(std::vector<double> matProperties)
    {
      m_volumeMaterialProperties = std::move(matProperties);
      m_volumeMaterial = Material(m_volumeMaterialProperties.at(0),
				  m_volumeMaterialProperties.at(1),
				  m_volumeMaterialProperties.at(2),
				  m_volumeMaterialProperties.at(3),
				  m_volumeMaterialProperties.at(4));
    }

    void setVolumeToBeamPipe(bool flag)
    {
      m_volumeToBeamPipe = flag;
    }
    
    void setLayerBuilder(std::shared_ptr<ILayerBuilder> layerBuilder)
    {
      m_layerBuilder = std::move(layerBuilder);
    }
    
    void setLayerEnvelopeR(double r)
    {
      m_layerEnvelopeR = r;
    }
    
    void setLayerEnvelopeZ(double z)
    {
      m_layerEnvelopeZ = z;
    }
    
    void setVolumeSignature(int sig)
    {
      m_volumeSignature = sig;
    }
    
  private:
    /** analyse the layer setup */        
    LayerSetup analyzeLayerSetup(const LayerVector* lVector) const;   

    std::shared_ptr<ITrackingVolumeHelper>     m_trackingVolumeHelper;       //!< the tracking volume creator for container volume creation

    std::string                            m_volumeName;                  //!< the name of the volume to be created
     
    std::vector< double >                  m_volumeDimension;             //!< The dimensions of the manually created world
    std::vector< double >                  m_volumeMaterialProperties;    //!< The material properties of the created world
    Material                               m_volumeMaterial;              //!< the world material  
    bool                                   m_volumeToBeamPipe;            //!< build the volume to the beam pipe

    std::shared_ptr<ILayerBuilder>              m_layerBuilder;                //!< needed to build layers within the volume
    double                                 m_layerEnvelopeR;              //!< the envelope covering the potential layers
    double                                 m_layerEnvelopeZ;              //!< the envelope covering the potential layers
      
    int                                    m_volumeSignature;             //!< the volume signature

  };

} // end of namespace

#endif // ACTS_GEOMETRYINTERFACES_WORLDVOLUMEBUILDER_H

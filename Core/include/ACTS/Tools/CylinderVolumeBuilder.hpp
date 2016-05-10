///////////////////////////////////////////////////////////////////
// CylinderVolumeBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYTOOLS_CYLINDERVOLUMEBUILDER_H
#define ACTS_GEOMETRYTOOLS_CYLINDERVOLUMEBUILDER_H 1

// Geometry module
#include "ACTS/Tools/ITrackingVolumeBuilder.hpp"
#include "ACTS/Material/Material.hpp"
#include "ACTS/Tools/ITrackingVolumeHelper.hpp"
#include "ACTS/Tools/ILayerBuilder.hpp"
#include "ACTS/Utilities/BinningType.hpp"

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
       
      @author Andreas.Salzburger@cern.ch   
  */
  
  class CylinderVolumeBuilder : public ITrackingVolumeBuilder{
    public:
      /** @struct Config 
          Configuration struct for this CylinderVolumeBuilder */
        struct Config {
            std::shared_ptr<ITrackingVolumeHelper>  trackingVolumeHelper;       //!< the tracking volume creator for container volume creation
            std::string                             volumeName;                  //!< the name of the volume to be created
            std::vector< double >                   volumeDimension;             //!< The dimensions of the manually created world
            Material                                volumeMaterial;              //!< the world material  
            bool                                    volumeToBeamPipe;            //!< build the volume to the beam pipe
            std::shared_ptr<ILayerBuilder>          layerBuilder;                //!< needed to build layers within the volume
            double                                  layerEnvelopeR;              //!< the envelope covering the potential layers
            double                                  layerEnvelopeZ;              //!< the envelope covering the potential layers
            int                                     volumeSignature;             //!< the volume signature 
        };   
        
        /** Constructor */
        CylinderVolumeBuilder(const Config& cvbConfig);
          
        /** Destructor */
        virtual ~CylinderVolumeBuilder();
          
        /** CylinderVolumeBuilder interface method - returns the volumes of Volumes */
        TrackingVolumePtr trackingVolume(TrackingVolumePtr insideVolume = nullptr,
        				                 VolumeBoundsPtr outsideBounds  = nullptr,
        				                 const LayerTriple* layerTriple  = nullptr,
        				                 const VolumeTriple* volumeTriple = nullptr) const override;
                               
       /** Set configuration method */
       void setConfiguration(const Config& cvbConfig);
      
       /** Get configuration method */
       Config getConfiguration() const;            
    
    protected:
        /** Configuration struct */
        Config m_config;                                   
      
    private:
        /** analyse the layer setup */        
        LayerSetup analyzeLayerSetup(const LayerVector lVector) const;   

  };
  
  /** Return the configuration object */    
  inline CylinderVolumeBuilder::Config CylinderVolumeBuilder::getConfiguration() const { return m_config; }

} // end of namespace

#endif // ACTS_GEOMETRYINTERFACES_WORLDVOLUMEBUILDER_H

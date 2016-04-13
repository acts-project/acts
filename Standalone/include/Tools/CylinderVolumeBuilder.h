///////////////////////////////////////////////////////////////////
// CylinderVolumeBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYTOOLS_CYLINDERVOLUMEBUILDER_H
#define ACTS_GEOMETRYTOOLS_CYLINDERVOLUMEBUILDER_H 1

// Core module
#include "CoreInterfaces/AlgToolBase.h"
// Geometry module
#include "GeometryInterfaces/ITrackingVolumeBuilder.h"
#include "Material/Material.h"
#include "GeometryInterfaces/ITrackingVolumeHelper.h"
#include "GeometryInterfaces/ILayerBuilder.h"
#include "GeometryUtils/BinningType.h"
// Gaudi
#include "GaudiKernel/ToolHandle.h"

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
  
    class CylinderVolumeBuilder : public AlgToolBase, virtual public ITrackingVolumeBuilder {

    public:
      /** Constructor */
      CylinderVolumeBuilder(const std::string&,const std::string&,const IInterface*);
      
      /** Destructor */
      virtual ~CylinderVolumeBuilder();
      
      /** AlgTool initialize method */
      virtual StatusCode initialize() override;
      
      /** AlgTool finalize method */
      virtual StatusCode finalize() override;
      
      /** CylinderVolumeBuilder interface method - returns the volumes of Volumes */
      TrackingVolumePtr trackingVolume(TrackingVolumePtr insideVolume = nullptr,
                                       VolumeBoundsPtr outsideBounds  = nullptr,
                                       const LayerTriple* layerTriple  = nullptr,
                                       const VolumeTriple* volumeTriple = nullptr) const override;

   private:
      /** analyse the layer setup */        
      LayerSetup analyzeLayerSetup(const LayerVector* lVector) const;   
     
      ToolHandle<ITrackingVolumeHelper>      m_trackingVolumeHelper;       //!< the tracking volume creator for container volume creation
     
      std::string                            m_volumeName;                  //!< the name of the volume to be created
     
      std::vector< double >                  m_volumeDimension;             //!< The dimensions of the manually created world
      std::vector< double >                  m_volumeMaterialProperties;    //!< The material properties of the created world
      Material                               m_volumeMaterial;              //!< the world material  
      bool                                   m_volumeToBeamPipe;            //!< build the volume to the beam pipe

      ToolHandle<ILayerBuilder>              m_layerBuilder;                //!< needed to build layers within the volume
      double                                 m_layerEnvelopeR;              //!< the envelope covering the potential layers
      double                                 m_layerEnvelopeZ;              //!< the envelope covering the potential layers
      
      int                                    m_volumeSignature;             //!< the volume signature

};

} // end of namespace

#endif // ACTS_GEOMETRYINTERFACES_WORLDVOLUMEBUILDER_H

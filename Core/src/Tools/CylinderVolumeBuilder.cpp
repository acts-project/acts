///////////////////////////////////////////////////////////////////
// CylinderVolumeBuilder.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Tools/CylinderVolumeBuilder.hpp"
#include "ACTS/Tools/ITrackingVolumeHelper.hpp"
#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Layers/CylinderLayer.hpp"
#include "ACTS/Layers/DiscLayer.hpp"
#include "ACTS/Volumes/CylinderVolumeBounds.hpp"
#include "ACTS/Surfaces/CylinderBounds.hpp"
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/MsgMacros.hpp"

// constructor
Acts::CylinderVolumeBuilder::CylinderVolumeBuilder(const Acts::CylinderVolumeBuilder::Config& cvbConfig):
  Acts::ITrackingVolumeBuilder(),
  m_config()
{
  setConfiguration(cvbConfig);
}    

// destructor
Acts::CylinderVolumeBuilder::~CylinderVolumeBuilder()
{}

// configuration
void Acts::CylinderVolumeBuilder::setConfiguration(const Acts::CylinderVolumeBuilder::Config& cvbConfig) {
  // @TODO check consistency    
  // copy the configuration 
  m_config = cvbConfig;
} 

std::shared_ptr<const Acts::TrackingVolume> Acts::CylinderVolumeBuilder::trackingVolume(TrackingVolumePtr insideVolume,
                                                                                        VolumeBoundsPtr outsideBounds,
                                                                                        const Acts::LayerTriple* layerTriple,
                                                                                        const VolumeTriple* volumeTriple) const
{
    // the return volume -----------------------------------------------------------------------------
    std::shared_ptr<const TrackingVolume> volume = nullptr;
    // used throughout
    TrackingVolumePtr nEndcap     = nullptr;
    TrackingVolumePtr barrel      = nullptr;
    TrackingVolumePtr pEndcap     = nullptr;    
        
    // get the full extend
    double volumeRmin  = 10e10;
    double volumeRmax  = -10e10;
    double volumeZmax  = 0.;
    
    // now analyize the layers that are provided -----------------------------------------------------
    LayerVector negativeLayers;
    LayerVector centralLayers;
    LayerVector positiveLayers;
    // - get the layers from a provided layer triple or from the layer builder
    if (layerTriple) {
        // the negative Layers
        negativeLayers = std::get<0>(*layerTriple);
        // the central Layers
        centralLayers  = std::get<1>(*layerTriple);
        // the positive Layer
        positiveLayers = std::get<2>(*layerTriple);
    }
    else {
        // the negative Layers
        negativeLayers = m_config.layerBuilder->negativeLayers();
        // the central Layers
        centralLayers  = m_config.layerBuilder->centralLayers();
        // the positive Layer
        positiveLayers = m_config.layerBuilder->positiveLayers();
    }
    // analyze the layers
    LayerSetup nLayerSetup = analyzeLayerSetup(negativeLayers);
    LayerSetup cLayerSetup = analyzeLayerSetup(centralLayers);
    LayerSetup pLayerSetup = analyzeLayerSetup(positiveLayers);
    // layer configuration --------------------------------------------------------------------------
    // dimensions
    double layerRmin = 10e10;
    double layerRmax = 0.;
    double layerZmax = 0.;
    // possbile configurations are:
    //
    // 111 - all     layers present 
    // 010 - central layers present
    //  0  - no layers present     
    int layerConfiguration = 0;
    if (nLayerSetup) {
        // negative layers are present 
        // MSG_DEBUG("Negative layers are present with r(min,max) / z(min,max) = " << nLayerSetup.rBoundaries << " / " << nLayerSetup.zBoundaries );
        takeSmaller(layerRmin,nLayerSetup.rBoundaries.first);
        takeBigger(layerRmax,nLayerSetup.rBoundaries.second);
        takeBigger(layerZmax,fabs(nLayerSetup.zBoundaries.first));
        // set the 100-digit for n present
        layerConfiguration += 100;
    }
    if (cLayerSetup) {
        // central layers are present
        // MSG_DEBUG("Central  layers are present with r(min,max) / z(min,max) = " << cLayerSetup.rBoundaries << " / " << cLayerSetup.zBoundaries << " layerRmin: " << layerRmin << " layerEnvelopeR: " << m_config.layerEnvelopeR);
        takeSmaller(layerRmin,cLayerSetup.rBoundaries.first);
        takeBigger(layerRmax,cLayerSetup.rBoundaries.second);
        takeBigger(layerZmax,fabs(cLayerSetup.zBoundaries.first));
        takeBigger(layerZmax,cLayerSetup.zBoundaries.second);
        // set the 10-digit for c present
        layerConfiguration += 10;
    }
    if (pLayerSetup) {
        // positive layers are present 
        // MSG_DEBUG("Positive layers are present with r(min,max) / z(min,max) = " << pLayerSetup.rBoundaries << " / " << pLayerSetup.zBoundaries << " layerRmin: " << layerRmin << " layerEnvelopeR: " << m_config.layerEnvelopeR);
        takeSmaller(layerRmin,pLayerSetup.rBoundaries.first);
        takeBigger(layerRmax,pLayerSetup.rBoundaries.second);
        takeBigger(layerZmax,pLayerSetup.zBoundaries.second);
        // set the 1-digit for p present
        layerConfiguration += 1;
    }
    
    if (layerConfiguration) {
        // MSG_DEBUG("Layer configuration estimated as " << layerConfiguration << " with r(min,max) / z(min,max) = " <<
    } else {
        // MSG_DEBUG("No layers present in this setup." );
    }
    
    // the inside volume dimensions ------------------------------------------------------------------
    double insideVolumeRmin = 0.;
    double insideVolumeRmax = 0.;
    double insideVolumeZmax = 0.;
    if (insideVolume){
       // cast to cylinder volume  
       const CylinderVolumeBounds* icvBounds = 
           dynamic_cast<const CylinderVolumeBounds*>(&(insideVolume->volumeBounds()));
       // cylindrical volume bounds are there
       if (icvBounds){ 
           // the outer radius of the inner volume
           insideVolumeRmin = icvBounds->innerRadius();
           insideVolumeRmax = icvBounds->outerRadius();
           insideVolumeZmax = insideVolume->center().z()+icvBounds->halflengthZ();
           // MSG_VERBOSE("Inner CylinderVolumeBounds provided from external builder, rMin/rMax/zMax = " << insideVolumeRmin << ", " << insideVolumeRmax << ", " << insideVolumeZmax);
       } else {
           // we need to bail out, the given volume is not cylindrical
           // MSG_ERROR("Given volume to wrap was not cylindrical. Bailing out.");
           // cleanup teh memory 
           negativeLayers.clear(); centralLayers.clear(); positiveLayers.clear();
           // return a null pointer, upstream builder will have to understand this
           return nullptr;
       }
    }
    
    // -------------------- outside boundary conditions -------------------------------------------------- 
    //check if we have outsideBounds 
    if (outsideBounds){
        const CylinderVolumeBounds* ocvBounds = dynamic_cast<const CylinderVolumeBounds*>(outsideBounds.get());
        // the cast to CylinderVolumeBounds needs to be successful
         if (ocvBounds){
             // get values from the out bounds 
             volumeRmin  = ocvBounds->innerRadius();
             volumeRmax  = ocvBounds->outerRadius();
             volumeZmax = ocvBounds->halflengthZ();
             // MSG_VERBOSE("Outer CylinderVolumeBounds provided from external builder, rMin/rMax/zMax = " << volumeRmin << ", " << volumeRmax << ", " << volumeZmax);
          } else {
             // MSG_ERROR("Non-cylindrical bounds given to the CylinderVolumeBuilder. Bailing out.");
             // cleanup teh memory 
             negativeLayers.clear(); centralLayers.clear(); positiveLayers.clear();
             // return a null pointer, upstream builder will have to understand this
             return nullptr;
         }
         // check if the outside bounds cover all the layers 
         if (layerConfiguration && (volumeRmin > layerRmin || volumeRmax < layerRmax || volumeZmax < layerZmax)){
             // MSG_ERROR("Given layer dimensions do not fit inside the provided volume bounds. Bailing out." << " volumeRmin: " << volumeRmin << " volumeRmax: " << volumeRmax << " layerRmin: " << layerRmin << " layerRmax: " << layerRmax << " volumeZmax: " << volumeZmax << " layerZmax: " << layerZmax);
             // cleanup teh memory 
             negativeLayers.clear(); centralLayers.clear(); positiveLayers.clear();
             // return a null pointer, upstream builder will have to understand this
             return nullptr; 
         }
    } else if (m_config.volumeDimension.size() > 2) {
        // cylinder volume
        // get values from the out bounds 
        volumeRmin  = m_config.volumeDimension[0];
        volumeRmax  = m_config.volumeDimension[1];
        volumeZmax  = m_config.volumeDimension[2];
        // MSG_VERBOSE("Outer CylinderVolumeBounds provided by configuration, rMin/rMax/zMax = " << volumeRmin << ", " << volumeRmax << ", " << volumeZmax);
        
    } else {
        // outside dimensions will have to be determined by the layer dimensions
        volumeRmin = m_config.volumeToBeamPipe ? 0. : layerRmin - m_config.layerEnvelopeR;
        volumeRmax = layerRmax + m_config.layerEnvelopeR;
        volumeZmax = layerZmax + m_config.layerEnvelopeZ;
        // from setup 
        // MSG_VERBOSE("Outer CylinderVolumeBounds estimated from layer setup, rMin/rMax/zMax = " << volumeRmin << ", " << volumeRmax << ", " << volumeZmax);
    }
    
    // -------------------- analyse the layer setups -------------------------------------------------- 
    TrackingVolumePtr negativeSector = nullptr;
    TrackingVolumePtr centralSector  = nullptr;
    TrackingVolumePtr positiveSector = nullptr;
    
    // wrapping condition
    // 0 - no wrapping
    // 1 - wrap central barrel plus endcap volumes around inside volume
    //   - (a) gap volumes will have to be created to extend to potential z extend (if required)
    // 2 - wrap full setup around inside volume (fitting)
    //   - (a) barrel without endcap volumes
    //   - (b) endcaps are present and their position in z is around the inside volume
    //   - (c) gap volumes will have to be created to extend to potential z extend (if required)
    int wrappingCondition = 0;
    //check if layers are present
    if (layerConfiguration){
        // screen output
        // MSG_DEBUG("Building Volume from layer configuration.");
        // barrel configuration
        double barrelRmin = 0.;
        double barrelRmax = 0.;
        double barrelZmax = 0.;
        // endcap configuration
        double endcapRmin = 0.;
        double endcapRmax = 0.;
        double endcapZmin = 0.;
        double endcapZmax = 0.;
        //if the containing volumes are given, get the boundaries of them
        if(volumeTriple) {
            VolumePtr nEndcapVolume = std::get<0>(*volumeTriple);
            VolumePtr barrelVolume = std::get<1>(*volumeTriple);
            VolumePtr endcapVolume = std::get<2>(*volumeTriple);
            if (barrelVolume) {
                const CylinderVolumeBounds* barrelBounds = dynamic_cast<const CylinderVolumeBounds*>(&(barrelVolume->volumeBounds()));
                barrelRmin = barrelBounds->innerRadius();
                barrelRmax = barrelBounds->outerRadius();
                barrelZmax = barrelVolume->center().z() + barrelBounds->halflengthZ();
                // MSG_VERBOSE("Outer Barrel bounds provided by configuration, rMin/rMax/zMax = " << barrelRmin << ", " << barrelRmax << ", " << barrelZmax);
            } else 
                // MSG_ERROR("No Barrel volume given for current hierarchy!");
            //check if end cap volumes are provided
            if (endcapVolume) {
                const CylinderVolumeBounds* endcapBounds = dynamic_cast<const CylinderVolumeBounds*>(&(endcapVolume->volumeBounds()));
                endcapRmin = endcapBounds->innerRadius();
                endcapRmax = endcapBounds->outerRadius();
                endcapZmin = fabs(endcapVolume->center().z())-endcapBounds->halflengthZ();
                endcapZmax = fabs(endcapVolume->center().z())+endcapBounds->halflengthZ();
                // MSG_VERBOSE("Outer Endcap bounds provided by configuration, rMin/rMax/zMin/zMax = " << endcapRmin << ", " << endcapRmax << ", " << endcapZmin << ", " << endcapZmax);
            }
            //now set the wrapping condition
            // wrapping condition can only be set if there's an inside volume
            if (insideVolume) {
                if (endcapVolume && endcapZmin < insideVolumeZmax) wrappingCondition = 1;
                else wrappingCondition = 2;
            }
        } else {
            //if no containing volumes are provided calculate the bounds from the layer configuration
            // wrapping condition can only be set if there's an inside volume
            if (insideVolume){
                if (insideVolumeRmax > volumeRmax || insideVolumeZmax > volumeZmax) {
                    // we need to bail out, the given volume does not fit around the other
                    // MSG_ERROR("Given layer dimensions do not fit around the provided inside volume. Bailing out." << "insideVolumeRmax: " << insideVolumeRmax << " layerRmin: " << layerRmin);
                    // cleanup teh memory
                    negativeLayers.clear(); centralLayers.clear(); positiveLayers.clear();
                    // return a null pointer, upstream builder will have to understand this
                    return nullptr;
                }
                if (pLayerSetup && pLayerSetup.zBoundaries.first < insideVolumeZmax) {
                    // set the barrel parameters
                    barrelRmin = insideVolumeRmin;
                    barrelRmax = volumeRmax;
                    barrelZmax = 0.5*(cLayerSetup.zBoundaries.second + pLayerSetup.zBoundaries.first);
                    // set the endcap parameters
                    endcapRmin = insideVolumeRmax;
                    endcapRmax = volumeRmax;
                    endcapZmin = barrelZmax;
                    endcapZmax = volumeZmax;
                    // set the wrapping condition
                    wrappingCondition = 1;
                }
                else {
                    // set the barrel parameters
                    barrelRmin = insideVolumeRmin;
                    barrelRmax = volumeRmax;
                    barrelZmax = cLayerSetup.zBoundaries.second < insideVolumeZmax ? insideVolumeZmax : volumeZmax;
                    // set the endcap parameters
                    endcapRmin = insideVolumeRmin;
                    endcapRmax = volumeRmax;
                    endcapZmin = barrelZmax;
                    endcapZmax = volumeZmax;
                    // set the wrapping condition
                    wrappingCondition = 2;
                }
            } else {
                // no inside volume is given, wrapping conditions remains 0
                barrelRmin = volumeRmin;
                barrelRmax = volumeRmax;
                barrelZmax = pLayerSetup ? 0.5*(cLayerSetup.zBoundaries.second+pLayerSetup.zBoundaries.first) : volumeZmax;
                // endcap parameters
                endcapRmin = volumeRmin;
                endcapRmax = volumeRmax;
                endcapZmin = barrelZmax;
                endcapZmax = volumeZmax;
            
            }
        }//else - no volume bounds given from translation
        
        // the barrel is created
        barrel = m_config.trackingVolumeHelper->createTrackingVolume(centralLayers,
                                                              m_config.volumeMaterial,
                                                              barrelRmin, barrelRmax,
                                                              -barrelZmax, barrelZmax,
                                                              m_config.volumeName+"::Barrel");
        
        
        // the negative endcap is created
        nEndcap = !negativeLayers.empty() ?
        m_config.trackingVolumeHelper->createTrackingVolume(negativeLayers,
                                                     m_config.volumeMaterial,
                                                     endcapRmin, endcapRmax,
                                                     -endcapZmax, -endcapZmin,
                                                     m_config.volumeName+"::NegativeEndcap") : nullptr;
        
        // the positive endcap is created
        pEndcap = !positiveLayers.empty() ?
        m_config.trackingVolumeHelper->createTrackingVolume(positiveLayers,
                                                     m_config.volumeMaterial,
                                                     endcapRmin, endcapRmax,
                                                     endcapZmin, endcapZmax,
                                                     m_config.volumeName+"::PositiveEndcap") : nullptr;
        
        
        // no wrapping condition
        if (wrappingCondition == 0){
            // we have endcap volumes
            if (nEndcap && pEndcap) {
                // a new barrel sector
                volume = m_config.trackingVolumeHelper->createContainerTrackingVolume({nEndcap, barrel, pEndcap});
            } else // just take the barrel as the return value
                volume = barrel;
            
        } else if (wrappingCondition == 1) {
            // a new barrel sector
            volume = m_config.trackingVolumeHelper->createContainerTrackingVolume({nEndcap, barrel, pEndcap});
            // now check if we need gaps as in 1
            if (fabs(insideVolumeZmax-volumeZmax) > 10e-5 ){
                // create the gap volumes
                // - negative side
                nEndcap = m_config.trackingVolumeHelper->createGapTrackingVolume(m_config.volumeMaterial,
                                                                          insideVolumeRmax, volumeRmax,
                                                                          -volumeZmax,-barrelZmax,
                                                                          1, false,
                                                                          m_config.volumeName+"::NegativeGap");
                // - positive side
                pEndcap = m_config.trackingVolumeHelper->createGapTrackingVolume(m_config.volumeMaterial,
                                                                          insideVolumeRmax, volumeRmax,
                                                                          barrelZmax, volumeZmax,
                                                                          1, false,
                                                                          m_config.volumeName+"::PositiveGap");
                // update the volume with the two sides
                insideVolume = m_config.trackingVolumeHelper->createContainerTrackingVolume({nEndcap, insideVolume, pEndcap});
            }
            // update the volume
            volume = m_config.trackingVolumeHelper->createContainerTrackingVolume({insideVolume, volume});
            
        } else if (wrappingCondition == 2){
            //create gap volumes if needed
            if (barrelZmax>insideVolumeZmax) {
                // create the gap volumes
                auto niGap = m_config.trackingVolumeHelper->createGapTrackingVolume(m_config.volumeMaterial,
                                                                             insideVolumeRmin, volumeRmin,
                                                                             -barrelZmax,-insideVolumeZmax,
                                                                             1, false,
                                                                             m_config.volumeName+"::InnerNegativeGap");
                
                auto piGap = m_config.trackingVolumeHelper->createGapTrackingVolume(m_config.volumeMaterial,
                                                                             insideVolumeRmin, volumeRmin,
                                                                             insideVolumeZmax, barrelZmax,
                                                                             1, false,
                                                                             m_config.volumeName+"::InnerPositiveGap");
                // pack into a new insideVolume
                insideVolume =  m_config.trackingVolumeHelper->createContainerTrackingVolume({niGap, insideVolume, piGap});
            }
            // create the container of the detector
            insideVolume = m_config.trackingVolumeHelper->createContainerTrackingVolume({insideVolume, barrel});
            volume = (nEndcap && pEndcap) ? m_config.trackingVolumeHelper->createContainerTrackingVolume({nEndcap, insideVolume, pEndcap}) : insideVolume;
        }

    } else if (outsideBounds){
        // screen output
        // MSG_DEBUG("Building Volume without layer configuration.");
        if (insideVolume && outsideBounds){
            // the barrel is created
            barrel = m_config.trackingVolumeHelper->createTrackingVolume({},
                                                                  m_config.volumeMaterial,
                                                                  insideVolumeRmin, volumeRmax,
                                                                  -insideVolumeZmax, insideVolumeZmax,
                                                                  m_config.volumeName+"::Barrel");
            // pack into the appropriate container
            volume = m_config.trackingVolumeHelper->createContainerTrackingVolume({insideVolume, barrel});
            // check if necap gaps are needed
            if (fabs(insideVolumeZmax-volumeZmax) > 10e-5){
                // the negative endcap is created
                nEndcap = m_config.trackingVolumeHelper->createTrackingVolume({},
                                                                       m_config.volumeMaterial,
                                                                       insideVolumeRmin, volumeRmax,
                                                                       -volumeZmax, -insideVolumeZmax,
                                                                       m_config.volumeName+"::NegativeEndcap");
                // the positive endcap is created
                pEndcap = m_config.trackingVolumeHelper->createTrackingVolume({},
                                                                       m_config.volumeMaterial,
                                                                       insideVolumeRmin, volumeRmax,
                                                                       insideVolumeZmax, volumeZmax,
                                                                       m_config.volumeName+"::PositiveEndcap");
                // pack into a the container
                volume = m_config.trackingVolumeHelper->createContainerTrackingVolume({nEndcap, barrel, pEndcap});
            }
            
        } else
            volume = TrackingVolume::create(nullptr, outsideBounds, m_config.volumeMaterial);
        
    } else {
      // MSG_ERROR("Neither layer configuration nor volume bounds given. Bailing out.");
    }
    
    // sign the volume
    volume->sign(GeometrySignature(m_config.volumeSignature));
    // now return what you have
    return volume;
}

Acts::LayerSetup Acts::CylinderVolumeBuilder::analyzeLayerSetup(const LayerVector lVector) const {
    // return object
    LayerSetup lSetup;
    // only if the vector is present it can actually be analyzed
    if (!lVector.empty()){
        // we have layers
        lSetup.present = true;
        for (auto& layer : lVector){
           // the thickness of the layer needs to be taken into account
           double thickness = layer->thickness();
           // get the center of the layer 
           const Vector3D& center = layer->surfaceRepresentation().center();
           // check if it is a cylinder layer
           const CylinderLayer* cLayer = dynamic_cast<const CylinderLayer*>(layer.get());
           if (cLayer){
               // set the binning to radial binning
               lSetup.binningValue = binR;
               // now we have access to all the information
               double rMinC  = cLayer->surfaceRepresentation().bounds().r()-0.5*thickness;
               double rMaxC  = cLayer->surfaceRepresentation().bounds().r()+0.5*thickness;
               double hZ     = cLayer->surfaceRepresentation().bounds().halflengthZ();
               takeSmaller(lSetup.rBoundaries.first,rMinC);
               takeBigger(lSetup.rBoundaries.second,rMaxC);
               takeSmaller(lSetup.zBoundaries.first,center.z()-hZ);
               takeBigger(lSetup.zBoundaries.second,center.z()+hZ);
           }
           // proceed further if it is a Disc layer
           const RadialBounds* dBounds = dynamic_cast<const RadialBounds*>(&(layer->surfaceRepresentation().bounds()));
           if (dBounds){
               // set the binning to radial binning
               lSetup.binningValue = binZ;
               // now we have access to all the information
               double rMinD =dBounds->rMin();
               double rMaxD =dBounds->rMax();
               double zMinD =  center.z()-0.5*thickness;
               double zMaxD =  center.z()+0.5*thickness;
               takeSmaller(lSetup.rBoundaries.first,rMinD);
               takeBigger(lSetup.rBoundaries.second,rMaxD);
               takeSmaller(lSetup.zBoundaries.first,zMinD);
               takeBigger(lSetup.zBoundaries.second,zMaxD);
               //!< @TODO check for Ring setup
           }
       }
    }    
    return lSetup;
} 

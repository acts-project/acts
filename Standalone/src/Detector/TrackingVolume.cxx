///////////////////////////////////////////////////////////////////
// TrackingVolume.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "Detector/TrackingVolume.h"
#include "Detector/GlueVolumesDescriptor.h"
#include "Detector/DetachedTrackingVolume.h"
#include "Volumes/VolumeBounds.h"
#include "Volumes/BoundaryCylinderSurface.h"
#include "Volumes/BoundaryDiscSurface.h"
#include "Volumes/BoundaryPlaneSurface.h"
#include "Volumes/BoundarySubtractedCylinderSurface.h"
#include "Volumes/BoundarySubtractedPlaneSurface.h"
#include "Volumes/CombinedVolumeBounds.h"
#include "Volumes/SubtractedVolumeBounds.h"
#include "Volumes/CylinderVolumeBounds.h"
#include "Volumes/SimplePolygonBrepVolumeBounds.h"
#include "Surfaces/Surface.h"
#include "Surfaces/CylinderSurface.h"
#include "Surfaces/DiscSurface.h"
#include "Surfaces/PlaneSurface.h"
#include "Surfaces/SubtractedCylinderSurface.h"
#include "Surfaces/SubtractedPlaneSurface.h"
#include "GeometryUtils/BinUtility.h"

// default constructor
Acts::TrackingVolume::TrackingVolume() :
  Volume(),
  Material(),
  m_motherVolume(),
  m_boundarySurfaces(nullptr),
  m_confinedLayers(nullptr),
  m_confinedVolumes(nullptr),
  m_confinedDetachedVolumes(nullptr),
  m_confinedDenseVolumes(nullptr),
  m_confinedArbitraryLayers(nullptr),
  m_glueVolumeDescriptor(nullptr),
  m_geometrySignature(Unsigned),
  m_geometryType(NumberOfGeometryTypes),
  m_name("undefined"),
  m_colorCode(20)
{}

// constructor for a container
Acts::TrackingVolume::TrackingVolume(std::shared_ptr<Transform3D> htrans,
                                    VolumeBoundsPtr volbounds,
                                    const std::shared_ptr<const TrackingVolumeArray> containedVolumeArray,
                                    const std::string& volumeName) :
  Volume(htrans, volbounds),
  Material(),
  m_motherVolume(nullptr),
  m_boundarySurfaces(nullptr),
  m_confinedLayers(nullptr),
  m_confinedVolumes(containedVolumeArray),
  m_confinedDetachedVolumes(nullptr),
  m_confinedDenseVolumes(nullptr),
  m_confinedArbitraryLayers(nullptr),
  m_glueVolumeDescriptor(nullptr),
  m_geometrySignature(Unsigned),
  m_geometryType(NumberOfGeometryTypes),
  m_name(volumeName),
  m_colorCode(20)
{
  createBoundarySurfaces();
  interlinkLayers();
}   

// constructor for arguments
Acts::TrackingVolume::TrackingVolume(std::shared_ptr<Transform3D> htrans,
                                    VolumeBoundsPtr volbounds,
                                    const Material& matprop,
                                    const LayerArray* staticLayerArray,
                                    const LayerVector* arbitraryLayerVector, 
                                    const TrackingVolumeArray* containedVolumeArray,
                                    const TrackingVolumeVector* denseVolumeVector,
                                    const DetachedVolumeVector* detachedVolumeVector,
                                    const std::string& volumeName) :
  Volume(htrans, volbounds),
  Material(matprop),
  m_motherVolume(nullptr),
  m_confinedLayers(staticLayerArray),
  m_confinedVolumes(containedVolumeArray),
  m_confinedDetachedVolumes(detachedVolumeVector),
  m_confinedDenseVolumes(denseVolumeVector),
  m_confinedArbitraryLayers(arbitraryLayerVector),
  m_glueVolumeDescriptor(nullptr),
  m_geometrySignature(Unsigned),
  m_geometryType(NumberOfGeometryTypes),
  m_name(volumeName),
  m_colorCode(20)
{
  createBoundarySurfaces();
  interlinkLayers();
}

Acts::TrackingVolume::TrackingVolume(const TrackingVolume& tvol,
                                    const Transform3D& shift, 
                                    const std::string& volumeName) :
  Volume(tvol, shift),
  Material(tvol),
  m_motherVolume(tvol.motherVolume()),
  m_confinedLayers(nullptr),
  m_confinedVolumes(nullptr),
  m_confinedDetachedVolumes(nullptr),
  m_confinedDenseVolumes(nullptr),
  m_confinedArbitraryLayers(nullptr),
  m_glueVolumeDescriptor(nullptr),
  m_geometrySignature(tvol.geometrySignature()),
  m_geometryType(tvol.geometryType()),
  m_name(volumeName),
  m_colorCode(20)
{    
    //< @TODO implement - requires cloneWithShift for BinUtility and an orderPosition() addon to GeometryObjects   
}                          

Acts::TrackingVolume::~TrackingVolume()
{
   delete m_boundarySurfaces;   
   delete m_confinedLayers;
   delete m_confinedDetachedVolumes;
   delete m_confinedDenseVolumes;
   delete m_confinedArbitraryLayers;
   delete m_glueVolumeDescriptor;
}

const Acts::Layer* Acts::TrackingVolume::associatedLayer(const Vector3D& gp) const
{  
  // confined static layers - highest hierarchy
  if (m_confinedLayers)
     return (m_confinedLayers->object(gp).get());
  
  // confined arbitrary   
  if (m_confinedArbitraryLayers) 
    for (auto& layer : (*m_confinedArbitraryLayers) ) 
      if ( layer->isOnLayer(gp) ) return layer.get();
  
  // return the null pointer 
  return nullptr;
}
                                                        
const Acts::TrackingVolume* Acts::TrackingVolume::associatedSubVolume(const Vector3D& gp) const
{
  // confined static volumes - highest hierarchy
  if (m_confinedVolumes) 
     return (m_confinedVolumes->object(gp).get());

  // if no static volumes are there, detached is next hierarchy
  if (m_confinedDetachedVolumes)
    for (auto& detachedVolume : (*m_confinedDetachedVolumes))  
        if (detachedVolume->trackingVolume()->inside(gp,0.001))
            return detachedVolume->trackingVolume();
      
  // if no static volumes or detached volumes are there, search for dense volumes
  if (m_confinedDenseVolumes)
    for (auto& denseVolume : (*m_confinedDenseVolumes)) 
        if (denseVolume->inside(gp,0.001))
            return denseVolume.get(); 
    
  // there is no lower sub structure
  return this;
}
      
const Acts::TrackingVolume* Acts::TrackingVolume::nextVolume(const Vector3D& gp, const Vector3D& dir, PropDirection pDir) const {
    // get the boundary surfaces & intersect them
    const TrackingVolume* nVolume = 0;
    // fix the direction once
    bool forceDir = ( pDir == alongMomentum || pDir == oppositeMomentum);
    double dirScalor = ( pDir == oppositeMomentum)  ? -1. : 1.;
    Vector3D cDir = dirScalor*dir;
    double pathLength = 10e10;
    // now loop through the and find the closest 
    auto bSurfaces = boundarySurfaces();
    for (auto& bSurfIter : bSurfaces ){
        // get the intersection soltuion
        Intersection sfI = bSurfIter->surfaceRepresentation().intersectionEstimate(gp, cDir, forceDir, true);
        if (sfI.valid && (sfI.pathLength*sfI.pathLength) < (pathLength*pathLength) ){
            // assign the next Volume
            PropDirection attachedDir = sfI.pathLength > 0. ? alongMomentum : oppositeMomentum;
            pathLength = sfI.pathLength;
            nVolume    = bSurfIter->attachedVolume(gp, cDir, attachedDir );
        }
    }
    return nVolume;
}
  
const Acts::DetachedVolumeVector* Acts::TrackingVolume::assocDetachedSubVolumes(const Vector3D& gp, double tol) const
{
  // create a new vector     
  DetachedVolumeVector* currVols = new DetachedVolumeVector;
  // get the volumes were the position is inside 
  if (m_confinedDetachedVolumes)
    for (auto& detachedVolume : (*m_confinedDetachedVolumes))
      if (detachedVolume->trackingVolume()->inside(gp,tol))
          currVols->push_back(detachedVolume);         
  // return the volumes that are inside 
  return currVols;
}

void Acts::TrackingVolume::addMaterial(const Material& mprop, float fact) const
{
  // assume the scaling factor refers to the volume scaling
  float flin = pow(fact,0.33); 
  //average X0
  double invX0 = X0>0. ? 1./X0 : 0.;
  double sum_invX0 = invX0 + flin/mprop.X0;
  X0 = 1./sum_invX0;
  //average L0
  double invL0 = L0>0. ? 1./L0 : 0.;
  double sum_invL0 = invL0 + flin/mprop.L0;
  L0 = 1./sum_invL0;
  //add density
  float rho1 = rho;
  rho += fact*mprop.rho; 
  // averageZ
  float n1 = Z>0. ? rho1/Z : 0.;  
  float n2 = fact*mprop.rho/mprop.Z;
  Z = rho/(n1+n2);   
  // averageA
  n1 = A>0. ? rho1/A : 0.;  
  n2 = fact*mprop.rho/mprop.A;
  A = rho/(n1+n2);   
  // zOverAtimesRho 
  zOaTr = Z/A*rho;  
  // mean energy loss (linear scaling)
  dEdX += flin*mprop.dEdX; 
}

void Acts::TrackingVolume::sign(GeometrySignature geosign, GeometryType geotype) const
{

  // never overwrite what is already signed, that's a crime
  if (m_geometrySignature == Unsigned) m_geometrySignature = geosign;
  m_geometryType = geotype;

  // confined static volumes
  if (m_confinedVolumes)
    for ( auto& volumesIter : (m_confinedVolumes->arrayObjects()) )
        volumesIter->sign(geosign, geotype);
  
  // same procedure for the detached volumes
  if (m_confinedDetachedVolumes)
     for (auto& volumesIter : (*m_confinedDetachedVolumes) )
        volumesIter->sign(geosign, geotype);
  
  // finally for confined dense volumes
  if (m_confinedDenseVolumes)
     for (auto& volumesIter : (*m_confinedDenseVolumes) )
       volumesIter->sign(geosign, geotype);
}

const std::vector< std::shared_ptr<const Acts::BoundarySurface<Acts::TrackingVolume> > >&  Acts::TrackingVolume::boundarySurfaces() const
{ return (*m_boundarySurfaces); }

void Acts::TrackingVolume::createBoundarySurfaces()
{
  // prepare the BoundarySurfaces
  m_boundarySurfaces = new std::vector< std::shared_ptr<const BoundarySurface<TrackingVolume> > >;
  
  // transform Surfaces To BoundarySurfaces
  const std::vector<const Surface*>* surfaces = Volume::volumeBounds().decomposeToSurfaces(m_transform);
  std::vector<const Surface*>::const_iterator surfIter = surfaces->begin();

  // counter to flip the inner/outer position for Cylinders
  unsigned int sfCounter = 0;
  unsigned int sfNumber  = surfaces->size();

  // memory optimisation
  m_boundarySurfaces->reserve(sfNumber+1);
  
  // identify Subtracted/CombinedVolumes
  const SubtractedVolumeBounds* subtrVol = dynamic_cast<const SubtractedVolumeBounds*> 
                                                                    (&( Volume::volumeBounds()));
  const CombinedVolumeBounds*   combVol  = dynamic_cast<const CombinedVolumeBounds*> 
                                                                    (&( Volume::volumeBounds()));
  bool subtr = (subtrVol) ? 1 : 0;
  bool comb  = (combVol)  ? 1 : 0;

  if (!subtr && !comb) {

    const SimplePolygonBrepVolumeBounds*   spbVol  = dynamic_cast<const SimplePolygonBrepVolumeBounds*> 
                                                                    (&( Volume::volumeBounds()));  

    for ( ; surfIter != surfaces->end(); ++surfIter){
      sfCounter++;

      TrackingVolume* in  = this;
      TrackingVolume* out = 0;
   
      // ST update: subtracted surfaces may appear in 'simple' volumes (SimplePolygonBrep...)
      const SubtractedPlaneSurface*      spsf = dynamic_cast<const SubtractedPlaneSurface*>(*surfIter);
      const PlaneSurface*      psf = dynamic_cast<const PlaneSurface*>(*surfIter);
      if (spsf) { 
        if (spbVol && sfCounter==1 ) {in = 0; out = this;}
	            m_boundarySurfaces->push_back(std::shared_ptr<const BoundarySurface<TrackingVolume> >
				      (new BoundarySubtractedPlaneSurface<TrackingVolume>(in, out, *spsf)));    
	    delete spsf; continue;
      } else if (psf){ m_boundarySurfaces->push_back(std::shared_ptr<const BoundarySurface<TrackingVolume> >
                                      (new BoundaryPlaneSurface<TrackingVolume>(in, out, *psf)));    
                     delete psf; continue;
      }        

      const DiscSurface*       dsf = dynamic_cast<const DiscSurface*>(*surfIter);
      if (dsf) {
          m_boundarySurfaces->push_back(std::shared_ptr<const BoundarySurface<TrackingVolume> >
              (new BoundaryDiscSurface<TrackingVolume>(in, out, *dsf)));
        delete dsf; continue;
      }

      const SubtractedCylinderSurface*  scsf = dynamic_cast<const SubtractedCylinderSurface*>(*surfIter);
      const CylinderSurface*   csf = dynamic_cast<const CylinderSurface*>(*surfIter);
      if (scsf) {
          TrackingVolume* inner = (sfCounter == 4 && sfNumber > 3) ? 0 : this;
          TrackingVolume* outer = (inner) ? 0 : this;
          m_boundarySurfaces->push_back(std::shared_ptr<const BoundarySurface< TrackingVolume> >
            (new BoundarySubtractedCylinderSurface<TrackingVolume>(inner, outer, *scsf)));
          delete scsf; continue;
      } else if (csf) {
          TrackingVolume* inner = (sfCounter == 4 && sfNumber > 3) ? 0 : this;
          TrackingVolume* outer = (inner) ? 0 : this;
          m_boundarySurfaces->push_back(std::shared_ptr<const BoundarySurface< TrackingVolume> >
            (new BoundaryCylinderSurface<TrackingVolume>(inner, outer, *csf)));
          delete csf; continue;
      } 
    }

  } else {

    const std::vector<bool> bOrient = subtrVol ? subtrVol->boundsOrientation() : combVol->boundsOrientation() ;
     
    for ( ; surfIter != surfaces->end(); ++surfIter){
      
      TrackingVolume* in  = bOrient[sfCounter] ? this : 0;
      TrackingVolume* out = bOrient[sfCounter] ? 0 : this;
      sfCounter++;
   
      const SubtractedPlaneSurface*      psf = dynamic_cast<const SubtractedPlaneSurface*>(*surfIter);
      if (psf) {
          m_boundarySurfaces->push_back(std::shared_ptr<const BoundarySurface<TrackingVolume> >
              (new BoundarySubtractedPlaneSurface<TrackingVolume>(in, out, *psf)));
          delete psf; continue;
      }        
      
      const SubtractedCylinderSurface*   csf = dynamic_cast<const SubtractedCylinderSurface*>(*surfIter);
      if (csf) {
          m_boundarySurfaces->push_back(std::shared_ptr<const BoundarySurface< TrackingVolume> >
              (new BoundarySubtractedCylinderSurface<TrackingVolume>(in, out, *csf)));
          delete csf; continue;
      }
    }
  }
  
  delete surfaces;
}

/** glue another tracking volume to this one */
void Acts::TrackingVolume::glueTrackingVolume(BoundarySurfaceFace bsfMine,
                                             std::shared_ptr<const TrackingVolume> neighbor,
                                             BoundarySurfaceFace bsfNeighbor) const 
{
    // find the connection of the two tracking volumes : binR returns the center except for cylindrical volumes
    Vector3D bPosition(binningPosition(binR));
    Vector3D distance = Vector3D(neighbor->binningPosition(binR)-bPosition);
    // glue to the face
    std::shared_ptr<const BoundarySurface<TrackingVolume> > bSurfaceMine = boundarySurfaces()[bsfMine]; 
    // @TODO - complex glueing could be possible with actual intersection for the normal vector 
    Vector3D normal = bSurfaceMine->surfaceRepresentation().normal(bPosition);
    // estimate the orientation
    BoundaryOrientation bOrientation = (normal.dot(distance) > 0.) ? outsideVolume : insideVolume;
    // the easy case : 
    // - no glue volume descriptors on either side 
    if (!m_glueVolumeDescriptor || !m_glueVolumeDescriptor->glueVolumes(bsfMine)){
        // the boundary orientation 
        bSurfaceMine->attachVolume(neighbor, bOrientation);
        // now set it to the neighbor volume - the optised way
        (*(neighbor->m_boundarySurfaces))[bsfNeighbor] =  bSurfaceMine; 
    }
}                                                

/** glue another tracking volume to this one */
void Acts::TrackingVolume::glueTrackingVolumes(BoundarySurfaceFace bsfMine,
                                              std::shared_ptr<const TrackingVolumeArray> neighbors,
                                              BoundarySurfaceFace bsfNeighbor) const 
{
    // find the connection of the two tracking volumes : binR returns the center except for cylindrical volumes    
    std::shared_ptr<const TrackingVolume> nRefVolume = neighbors->arrayObjects()[0];
    // get the distance 
    Vector3D bPosition(binningPosition(binR));
    Vector3D distance = Vector3D(nRefVolume->binningPosition(binR)-bPosition);
    // take the normal at the binning positio
    std::shared_ptr<const BoundarySurface<TrackingVolume> > bSurfaceMine = boundarySurfaces()[bsfMine];
    // @TODO - complex glueing could be possible with actual intersection for the normal vector 
    Vector3D normal = bSurfaceMine->surfaceRepresentation().normal(bPosition);
    // estimate the orientation
    BoundaryOrientation bOrientation = (normal.dot(distance) > 0.) ? outsideVolume : insideVolume;
    // the easy case : 
    // - no glue volume descriptors on either side 
    if (!m_glueVolumeDescriptor || !m_glueVolumeDescriptor->glueVolumes(bsfMine)){
        // the boundary orientation 
        bSurfaceMine->attachVolumeArray(neighbors, bOrientation);
        // now set it to the neighbor volumes - the optised way
        for (auto& nVolume : neighbors->arrayObjects())        
            (*(nVolume->m_boundarySurfaces))[bsfNeighbor] =  bSurfaceMine; 
    }
} 

/** update the boundary surface after glueing */
void Acts::TrackingVolume::updateBoundarySurface(BoundarySurfaceFace bsf, std::shared_ptr<const BoundarySurface<TrackingVolume> > bs) const
{
    (*m_boundarySurfaces)[bsf] = bs;
}    
    
void Acts::TrackingVolume::registerGlueVolumeDescriptor(GlueVolumesDescriptor* gvd) const 
{ 
    delete m_glueVolumeDescriptor;
    m_glueVolumeDescriptor = gvd;
}

const Acts::GlueVolumesDescriptor& Acts::TrackingVolume::glueVolumesDescriptor() const 
{
    if (!m_glueVolumeDescriptor)
        m_glueVolumeDescriptor = new GlueVolumesDescriptor;
    return (*m_glueVolumeDescriptor);
}

void Acts::TrackingVolume::synchronizeLayers(double envelope) const {

  // case a : Layers exist
  // msgstream << MSG::VERBOSE << "  -> synchronizing Layer dimensions of TrackingVolume '" << volumeName() << "'." << endreq;     
    
  if (m_confinedLayers){
    // msgstream << MSG::VERBOSE << "  ---> working on " << m_confinedLayers->arrayObjects().size() << " (material+navigation) layers." << endreq;
    for (auto& clayIter : m_confinedLayers->arrayObjects())
        if (clayIter){
          // @TODO implement syncrhonize layer     
          //  if (clayIter->surfaceRepresentation().type() == Surface::Cylinder && !(center().isApprox(clayIter->surfaceRepresentation().center())) )
          //      clayIter->resizeAndRepositionLayer(volumeBounds(),center(),envelope);
          //  else 
          //      clayIter->resizeLayer(volumeBounds(),envelope);
        }  // else
            // msgstream << MSG::WARNING << "  ---> found 0 pointer to layer, indicates problem." << endreq;
  }

  // case b : container volume -> step down
  if (m_confinedVolumes){
    // msgstream << MSG::VERBOSE << "  ---> no confined layers, working on " << m_confinedVolumes->arrayObjects().size() << " confined volumes." << endreq;
    for (auto& cVolumesIter : m_confinedVolumes->arrayObjects())
        cVolumesIter->synchronizeLayers(envelope);
  }
   
}

void Acts::TrackingVolume::interlinkLayers() {
    
   // std::cout << "TrackingVolume::interlinkLayers0" << std::endl;
  if (m_confinedLayers){
   //   std::cout << "TrackingVolume::interlinkLayers1" << std::endl;
    auto& layers = m_confinedLayers->arrayObjects();
  //  std::cout << "TrackingVolume::interlinkLayers2" << std::endl;
    // forward register the last one as the previous one
    //  first <- | -> second, first <- | -> second, first <- | -> second
    const Layer* lastLayer = nullptr;
    for (auto& layerPtr : layers)
      {
   //       std::cout << "TrackingVolume::interlinkLayers03" << std::endl;
        // register the layers
        (*layerPtr).m_nextLayerUtility = m_confinedLayers->binUtility();
   //       std::cout << "TrackingVolume::interlinkLayers4" << std::endl;
        (*layerPtr).m_nextLayers.first = lastLayer;
   //       std::cout << "TrackingVolume::interlinkLayers05" << std::endl;
        // register the volume
        (*layerPtr).encloseTrackingVolume(*this);
   //       std::cout << "TrackingVolume::interlinkLayers06" << std::endl;
        // remember the last layer
        lastLayer = layerPtr.get();
   //       std::cout << "TrackingVolume::interlinkLayers07" << std::endl;
      }
    // backward loop
    lastLayer = nullptr;
    for ( auto layerIter = layers.rbegin(); layerIter != layers.rend(); ++layerIter )
      {
   //       std::cout << "TrackingVolume::interlinkLayers8" << std::endl;
        // set the other next volume
        (**layerIter).m_nextLayers.second = lastLayer;
   //       std::cout << "TrackingVolume::interlinkLayers9" << std::endl;
        lastLayer = (*layerIter).get();
   //       std::cout << "TrackingVolume::interlinkLayers10" << std::endl;
      }
  }
}

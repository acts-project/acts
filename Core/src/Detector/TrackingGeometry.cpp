///////////////////////////////////////////////////////////////////
// TrackingGeometry.cpp, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Detector/DetachedTrackingVolume.hpp"
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Surfaces/PerigeeSurface.hpp"

Acts::TrackingGeometry::TrackingGeometry(TrackingVolumePtr highestVolume)
: m_world(highestVolume),
m_beam(std::make_unique<const Acts::PerigeeSurface>())
{
    // register all the TrackingVolumes 
    if (m_world) 
        registerTrackingVolumes(*m_world.get());
}

Acts::TrackingGeometry::~TrackingGeometry()
{}

const Acts::TrackingVolume* Acts::TrackingGeometry::lowestTrackingVolume(const Acts::Vector3D& gp) const
{
    const Acts::TrackingVolume* searchVolume  = m_world.get();
    const Acts::TrackingVolume* currentVolume = nullptr;
    while (currentVolume != searchVolume && searchVolume) {
        currentVolume = searchVolume;
        searchVolume  = searchVolume->associatedSubVolume(gp);
    }
    return currentVolume;
}

const Acts::DetachedVolumeVector* Acts::TrackingGeometry::lowestDetachedTrackingVolumes(const Acts::Vector3D& gp) const
{
    double tol = 0.001;
    const Acts::TrackingVolume* currentVolume = lowestStaticTrackingVolume(gp);
    if (currentVolume) return currentVolume->assocDetachedSubVolumes(gp,tol);
    return nullptr;
}

const Acts::TrackingVolume* Acts::TrackingGeometry::lowestStaticTrackingVolume(const Acts::Vector3D& gp) const
{
    const Acts::TrackingVolume* searchVolume  = m_world.get();
    const Acts::TrackingVolume* currentVolume = nullptr;
    while (currentVolume != searchVolume && searchVolume ) {
        currentVolume = searchVolume;
        if (!(searchVolume->confinedDetachedVolumes()) ) 
            searchVolume  = searchVolume->associatedSubVolume(gp);
    }
    return currentVolume;
}

void Acts::TrackingGeometry::registerTrackingVolumes(const Acts::TrackingVolume& tvol, const Acts::TrackingVolume* mvol, int lvl)
{
    int sublvl = lvl+1;
    std::string indent = "";
    for (int l=0; l<lvl; ++l, indent += "  ");

    tvol.setMotherVolume(mvol);
    
    m_trackingVolumes[tvol.volumeName()] = (&tvol);
    const Acts::TrackingVolumeArray* confinedVolumes = tvol.confinedVolumes();
    if (confinedVolumes){
        for (auto& volumesIter: confinedVolumes->arrayObjects())
            if (volumesIter) registerTrackingVolumes(*volumesIter, &tvol, sublvl);
    }

    const Acts::TrackingVolumeVector* confinedDenseVolumes = tvol.confinedDenseVolumes();
    if (confinedDenseVolumes){
        for (auto& volumesIter : (*confinedDenseVolumes) )
            if (volumesIter) registerTrackingVolumes(*volumesIter, &tvol, sublvl);
    }
    
    /** should detached tracking volumes be part of the tracking geometry ? */
    const Acts::DetachedVolumeVector* confinedDetachedVolumes = tvol.confinedDetachedVolumes();
    if (confinedDetachedVolumes){
        for (auto& volumesIter : (*confinedDetachedVolumes) )
            if (volumesIter && tvol.inside(volumesIter->trackingVolume()->center(),0.) )
                registerTrackingVolumes(*(volumesIter->trackingVolume()), &tvol, sublvl);
    }
    
}

//@TODO change to BoundaryCheck
bool Acts::TrackingGeometry::atVolumeBoundary( const Acts::Vector3D& gp, const Acts::TrackingVolume* vol, double) const
{
    bool isAtBoundary = false;
    if (!vol) return isAtBoundary;
    for (auto& bSurface: vol->boundarySurfaces()) {
        const Acts::Surface& surf = bSurface->surfaceRepresentation();
        if ( surf.isOnSurface(gp,true)  ) isAtBoundary = true;
    }
    return isAtBoundary;
} 

/** check position at volume boundary + navigation link */
//@TODO change to BoundaryCheck
bool Acts::TrackingGeometry::atVolumeBoundary(const Acts::Vector3D& gp, const Acts::Vector3D& mom, const TrackingVolume* vol, 
        const TrackingVolume*& nextVol, Acts::PropDirection dir, double) const
{
    bool isAtBoundary = false;
    nextVol = 0;
    if (!vol) return isAtBoundary;
    for (auto& bSurface: vol->boundarySurfaces()) {
        const Acts::Surface& surf = bSurface->surfaceRepresentation();
        if ( surf.isOnSurface(gp,true) ) {
            isAtBoundary = true;
            const Acts::TrackingVolume* attachedVol = bSurface->attachedVolume(gp,mom,dir);
            if (!nextVol && attachedVol) nextVol=attachedVol;
        }
    }
    return isAtBoundary;
} 

const Acts::TrackingVolume* Acts::TrackingGeometry::highestTrackingVolume() const
{ return (m_world.get()); }

void Acts::TrackingGeometry::sign(GeometrySignature geosit, GeometryType geotype) const
{ m_world->sign(geosit, geotype); }

const Acts::TrackingVolume* Acts::TrackingGeometry::trackingVolume(const std::string& name) const
{
    std::map<const std::string, const TrackingVolume*>::const_iterator sVol = m_trackingVolumes.begin();
    sVol = m_trackingVolumes.find(name);
    if (sVol != m_trackingVolumes.end()) return(sVol->second);
    return nullptr;
}

const Acts::Layer* Acts::TrackingGeometry::associatedLayer(const Acts::Vector3D& gp) const
{
    const TrackingVolume* lowestVol = (lowestTrackingVolume(gp));
    return lowestVol->associatedLayer(gp);
}

void Acts::TrackingGeometry::registerBeamTube(std::unique_ptr<const Acts::PerigeeSurface> beam) const
{
    m_beam = std::move(beam);
}






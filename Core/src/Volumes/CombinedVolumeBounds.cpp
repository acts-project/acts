///////////////////////////////////////////////////////////////////
// CombinedVolumeBounds.cpp, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "ACTS/Volumes/CombinedVolumeBounds.hpp"
#include "ACTS/Volumes/SubtractedVolumeBounds.hpp"
#include "ACTS/Volumes/VolumeExcluder.hpp"
#include "ACTS/Volumes/CylinderVolumeBounds.hpp"
#include "ACTS/Volumes/SimplePolygonBrepVolumeBounds.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/SubtractedPlaneSurface.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Surfaces/SubtractedCylinderSurface.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "ACTS/Surfaces/DiscSurface.hpp"
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Surfaces/EllipseBounds.hpp"
// STD/STL
#include <iostream>
#include <cmath>
#include <stdexcept>

Acts::CombinedVolumeBounds::CombinedVolumeBounds() :
 VolumeBounds(),
 m_first(nullptr),
 m_second(nullptr),
 m_intersection(false),
 m_boundsOrientation()
{}

Acts::CombinedVolumeBounds::CombinedVolumeBounds(Volume* vol1, Volume* vol2, bool intersection ) :
 VolumeBounds(),
 m_first(vol1),
 m_second(vol2),
 m_intersection(intersection),
 m_boundsOrientation()
{}

Acts::CombinedVolumeBounds::CombinedVolumeBounds(const Acts::CombinedVolumeBounds& bobo) :
 VolumeBounds(),
 m_first(bobo.m_first),
 m_second(bobo.m_second),
 m_intersection(bobo.m_intersection),
 m_boundsOrientation()
{
  m_boundsOrientation.resize(bobo.m_boundsOrientation.size());
  for (unsigned int i=0;i<bobo.m_boundsOrientation.size();i++) m_boundsOrientation[i]=bobo.m_boundsOrientation[i];
}

Acts::CombinedVolumeBounds::~CombinedVolumeBounds()
{
  m_boundsOrientation.clear();
  delete m_first;
  delete m_second;
}

Acts::CombinedVolumeBounds& Acts::CombinedVolumeBounds::operator=(const Acts::CombinedVolumeBounds& bobo)
{
  if (this!=&bobo){
    m_first          = bobo.m_first;
    m_second         = bobo.m_second;
    m_intersection   = bobo.m_intersection;
    m_boundsOrientation = bobo.m_boundsOrientation;
    m_boundsOrientation.resize(bobo.m_boundsOrientation.size());
    for (unsigned int i=0;i<bobo.m_boundsOrientation.size();i++) m_boundsOrientation[i]=bobo.m_boundsOrientation[i];
 }
  return *this;
}

const std::vector<const Acts::Surface*>* Acts::CombinedVolumeBounds::decomposeToSurfaces(std::shared_ptr<Acts::Transform3D> transformPtr) const
{

    std::vector<const Acts::Surface*>* retsf = new std::vector<const Acts::Surface*>;

    const Acts::CylinderVolumeBounds*   cylVol = dynamic_cast<const Acts::CylinderVolumeBounds*> (&(m_first->volumeBounds()));
    const Acts::SimplePolygonBrepVolumeBounds*   spbVol = dynamic_cast<const Acts::SimplePolygonBrepVolumeBounds*> (&(m_first->volumeBounds()));
    const Acts::CombinedVolumeBounds*   comVol = dynamic_cast<const Acts::CombinedVolumeBounds*> (&(m_first->volumeBounds()));
    const Acts::SubtractedVolumeBounds* subVol = dynamic_cast<const Acts::SubtractedVolumeBounds*> (&(m_first->volumeBounds()));

    // get surfaces for first boundaries
    Acts::Transform3D transf = ( transformPtr == nullptr) ? Acts::Transform3D::Identity() : (*transformPtr.get());

    std::shared_ptr<Acts::Transform3D> firstTransform(new Acts::Transform3D(transf*m_first->transform()));
    const std::vector<const Acts::Surface*>* firstSurfaces = m_first->volumeBounds().decomposeToSurfaces(firstTransform);

    // get surfaces for second boundaries
    std::shared_ptr<Acts::Transform3D> secondTransfrom(new Acts::Transform3D(transf*m_second->transform()));
    const std::vector<const Acts::Surface*>* secondSurfaces = m_second->volumeBounds().decomposeToSurfaces(secondTransfrom);
    unsigned int nSurf = firstSurfaces->size() + secondSurfaces->size();
    m_boundsOrientation.resize(nSurf);


    std::vector<unsigned int> subtrSecond;

    // loop over surfaces; convert disc surface to a plane surface using elliptic bounds
    for (unsigned int out=0; out < firstSurfaces->size(); out++) {
        //
        const SubtractedPlaneSurface* splo = dynamic_cast<const SubtractedPlaneSurface*> ((*firstSurfaces)[out]);
        const PlaneSurface* plo = dynamic_cast<const PlaneSurface*> ((*firstSurfaces)[out]);
        const SubtractedCylinderSurface* sclo = dynamic_cast<const SubtractedCylinderSurface*> ((*firstSurfaces)[out]);
        const CylinderSurface* clo = dynamic_cast<const CylinderSurface*> ((*firstSurfaces)[out]);
        const DiscSurface* dlo = dynamic_cast<const DiscSurface*> ((*firstSurfaces)[out]);

        // resolve bounds orientation : copy from combined/subtracted, swap inner cyl, swap bottom spb
        if (comVol) m_boundsOrientation[out]=comVol->boundsOrientation()[out];
        else if (subVol) m_boundsOrientation[out]=subVol->boundsOrientation()[out];
        else if (cylVol && clo && out==3 ) m_boundsOrientation[out] = false;
        else if (spbVol && out==0 ) m_boundsOrientation[out] = false;
        else m_boundsOrientation[out] = true;

        Acts::Volume* secondSub = createSubtractedVolume((*firstSurfaces)[out]->transform().inverse()*transf, m_second);

        if ( sclo || splo ) {
            bool shared = false;
            std::shared_ptr<Acts::AreaExcluder> vEx;
            if (splo) {
                vEx = splo->subtractedVolume();
                shared   = splo->shared();
            }
            if (sclo) {
                vEx = sclo->subtractedVolume();
                shared   = sclo->shared();
            }
            //vEx.addRef();
            const Acts::VolumeExcluder* volExcl = dynamic_cast<const Acts::VolumeExcluder*> (vEx.get());
            if (!volExcl) throw std::logic_error("Not a VolumeExcluder");

            Acts::Volume* firstSub = new Acts::Volume(*volExcl->volume());

            Acts::Volume* comb_sub = 0;
            if (!shared && !m_intersection) comb_sub = new Acts::Volume(nullptr,new Acts::CombinedVolumeBounds(secondSub,firstSub,m_intersection));
            if (!shared &&  m_intersection) comb_sub = new Acts::Volume(nullptr,new Acts::SubtractedVolumeBounds(secondSub,firstSub));
            if ( shared &&  m_intersection) comb_sub = new Acts::Volume(nullptr,new Acts::CombinedVolumeBounds(secondSub,firstSub,m_intersection));
            if ( shared && !m_intersection) comb_sub = new Acts::Volume(nullptr,new Acts::SubtractedVolumeBounds(firstSub,secondSub));
            Acts::VolumeExcluder* volEx = new Acts::VolumeExcluder(comb_sub);
            bool new_shared = shared;
            if (m_intersection) new_shared = true;
            if (splo) retsf->push_back(new Acts::SubtractedPlaneSurface(*splo,volEx,new_shared));
            if (sclo) retsf->push_back(new Acts::SubtractedCylinderSurface(*sclo,volEx,new_shared));

        }
        else if (plo || clo || dlo) {
            Acts::VolumeExcluder* volEx = new Acts::VolumeExcluder(secondSub);
            if (plo) retsf->push_back(new Acts::SubtractedPlaneSurface(*plo,volEx,m_intersection));
            if (clo) retsf->push_back(new Acts::SubtractedCylinderSurface(*clo,volEx,m_intersection));
            if (dlo) {
                const RadialBounds* db = dynamic_cast<const RadialBounds*> (&(dlo->bounds()));
                if (!db) throw std::logic_error("Not RadialBounds");

                EllipseBounds* eb = new EllipseBounds(db->rMin(),db->rMin(),db->rMax(),db->rMax(),db->halfPhiSector());
                plo = new PlaneSurface(std::shared_ptr<Acts::Transform3D>(new Acts::Transform3D(dlo->transform())),eb);
                retsf->push_back(new Acts::SubtractedPlaneSurface(*plo,volEx,m_intersection));
                delete plo;
            }
        }
        else {
            throw std::runtime_error("Unhandled surface in CombinedVolumeBounds::decomposeToSurfaces.");
        }
    }

    cylVol = dynamic_cast<const Acts::CylinderVolumeBounds*> (&(m_second->volumeBounds()));
    spbVol = dynamic_cast<const Acts::SimplePolygonBrepVolumeBounds*> (&(m_second->volumeBounds()));
    comVol = dynamic_cast<const Acts::CombinedVolumeBounds*> (&(m_second->volumeBounds()));
    subVol = dynamic_cast<const Acts::SubtractedVolumeBounds*> (&(m_second->volumeBounds()));
    unsigned int nOut = firstSurfaces->size();

    for (unsigned int in=0; in< secondSurfaces->size(); in++) {
        //
        const SubtractedPlaneSurface* spli = dynamic_cast<const SubtractedPlaneSurface*> ((*secondSurfaces)[in]);
        const PlaneSurface* pli = dynamic_cast<const PlaneSurface*> ((*secondSurfaces)[in]);
        const SubtractedCylinderSurface* scli = dynamic_cast<const SubtractedCylinderSurface*> ((*secondSurfaces)[in]);
        const CylinderSurface* cli = dynamic_cast<const CylinderSurface*> ((*secondSurfaces)[in]);
        const DiscSurface* dli = dynamic_cast<const DiscSurface*> ((*secondSurfaces)[in]);

        // resolve bounds orientation : copy from combined/subtracted, swap inner cyl, swap bottom spb
        if (comVol) m_boundsOrientation[nOut+in]=comVol->boundsOrientation()[in];
        else if (subVol) m_boundsOrientation[nOut+in]=subVol->boundsOrientation()[in];
        else if (cylVol && cli && in==3 ) m_boundsOrientation[nOut+in]= false;
        else if (spbVol && in==0 ) m_boundsOrientation[nOut+in]= false;
        else m_boundsOrientation[nOut+in]= true;

        Acts::Volume* firstSub = createSubtractedVolume((*secondSurfaces)[in]->transform().inverse()*transf, m_first);
        if ( scli || spli ) {
            bool shared = false;
            std::shared_ptr<Acts::AreaExcluder> vEx;
            if (spli) {
                vEx = spli->subtractedVolume();
                shared   = spli->shared();
            }
            if (scli) {
                vEx = scli->subtractedVolume();
                shared   = scli->shared();
            }
            //vEx.addRef();
            const Acts::VolumeExcluder* volExcl = dynamic_cast<const Acts::VolumeExcluder*> (vEx.get());
            if (!volExcl) throw std::logic_error("Not a VolumeExcluder");
            Acts::Volume* secondSub = new Acts::Volume(*volExcl->volume());

            Acts::Volume* comb_sub = 0;
            if (!shared && !m_intersection) comb_sub = new Acts::Volume(nullptr,new Acts::CombinedVolumeBounds(firstSub,secondSub,m_intersection));
            if (!shared &&  m_intersection) comb_sub = new Acts::Volume(nullptr,new Acts::SubtractedVolumeBounds(firstSub,secondSub));
            if ( shared &&  m_intersection) comb_sub = new Acts::Volume(nullptr,new Acts::CombinedVolumeBounds(firstSub,secondSub,m_intersection));
            if ( shared && !m_intersection) comb_sub = new Acts::Volume(nullptr,new Acts::SubtractedVolumeBounds(secondSub,firstSub));
            Acts::VolumeExcluder* volEx = new Acts::VolumeExcluder(comb_sub);
            bool new_shared = shared;
            if (m_intersection) new_shared = true;
            if (spli) retsf->push_back(new Acts::SubtractedPlaneSurface(*spli,volEx,new_shared));
            if (scli) retsf->push_back(new Acts::SubtractedCylinderSurface(*scli,volEx,new_shared));

        }
        else if (pli || cli || dli) {
            Acts::VolumeExcluder* volEx = new Acts::VolumeExcluder(firstSub);
            if (pli) retsf->push_back(new Acts::SubtractedPlaneSurface(*pli,volEx,m_intersection));
            if (cli) retsf->push_back(new Acts::SubtractedCylinderSurface(*cli,volEx,m_intersection));
            if (dli) {
                const RadialBounds* db = dynamic_cast<const RadialBounds*> (&(dli->bounds()));
                if (!db) throw std::logic_error("Not RadialBounds");

                EllipseBounds* eb = new EllipseBounds(db->rMin(),db->rMin(),db->rMax(),db->rMax(),db->halfPhiSector());
                pli = new PlaneSurface(std::shared_ptr<Acts::Transform3D>(new Acts::Transform3D(dli->transform())),eb);
                retsf->push_back(new Acts::SubtractedPlaneSurface(*pli,volEx,m_intersection));
                delete pli;
            }
        }
        else {
            throw std::runtime_error("Unhandled surface in CombinedVolumeBounds::decomposeToSurfaces.");
        }
    }

    for (size_t i=0; i < firstSurfaces->size(); i++)
        delete (*firstSurfaces)[i];
    for (size_t i=0; i < secondSurfaces->size(); i++)
        delete (*secondSurfaces)[i];
    delete firstSurfaces;
    delete secondSurfaces;

    return retsf;

}

// ostream operator overload
std::ostream& Acts::CombinedVolumeBounds::dump( std::ostream& sl ) const
{
    std::stringstream temp_sl;
    temp_sl << std::setiosflags(std::ios::fixed);
    temp_sl << std::setprecision(7);
    temp_sl << "Acts::CombinedVolumeBounds: first,second ";
    sl << temp_sl.str();
    m_first->volumeBounds().dump(sl);
    m_second->volumeBounds().dump(sl);
    return sl;
}

Acts::Volume* Acts::CombinedVolumeBounds::createSubtractedVolume(const Acts::Transform3D& transf, Acts::Volume* subtrVol) const
{
  Acts::Volume* subVol = 0;
  if (!subtrVol) return subVol;

  subVol = new Acts::Volume( *subtrVol, transf );

  return subVol;
}






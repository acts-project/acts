///////////////////////////////////////////////////////////////////
// SimplePolygonBrepVolumeBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_VOLUMES_SIMPLEPOLYGONBREPVOLUMEBOUNDS_H
#define ACTS_VOLUMES_SIMPLEPOLYGONBREPVOLUMEBOUNDS_H 1

#include "ACTS/Utilities/Definitions.h"
#include "ACTS/Volumes/VolumeBounds.h"

namespace Acts {

  class Surface;
  class PlaneSurface;
  class RectangleBounds;


  /**
   @class SimplePolygonBrepVolumeBounds

   Bounds for the exact transcript of the GeoSimplePolygonBrep; volume defined by combination of symm.trapezoids

    BoundarySurfaceFace [index]:

        - negativeFaceXY     [0] : Acts::SubtractedPlaneSurface,
                                   parallel to \f$ xy \f$ plane at negative \f$ z \f$
        - positiveFaceXY     [1] : Acts::SubtractedPlaneSurface,
                                   parallel to \f$ xy \f$ plane at positive \f$ z \f$
        - face [2... n+1] : Rectangular  Acts::PlaneSurface

    @author sarka.todorova@cern.ch , marcin.wolter@cern.ch
    */

 class SimplePolygonBrepVolumeBounds : public VolumeBounds {

  public:
    /**Default Constructor*/
    SimplePolygonBrepVolumeBounds();

    /**Constructor - generic case (from float)*/
    SimplePolygonBrepVolumeBounds(std::vector<std::pair<float,float> > xyvtx, float hlengthz);

    /**Constructor - generic case (from double)*/
    SimplePolygonBrepVolumeBounds(std::vector<std::pair<double,double> > xyvtx, double hlengthz);

    /**Copy Constructor */
    SimplePolygonBrepVolumeBounds(const SimplePolygonBrepVolumeBounds& bobo);

    /**Destructor */
    virtual ~SimplePolygonBrepVolumeBounds();

    /**Assignment operator*/
    SimplePolygonBrepVolumeBounds& operator=(const SimplePolygonBrepVolumeBounds& bobo);

    /**Virtual constructor */
    SimplePolygonBrepVolumeBounds* clone() const override;

    /**This method checks if position in the 3D volume frame is inside the volume*/
    bool inside(const Vector3D& , double tol=0.) const override;

    /** Method to decompose the Bounds into Surfaces */
    const std::vector<const Acts::Surface*>* decomposeToSurfaces(std::shared_ptr<Transform3D> transformPtr) const override;

    /**This method returns the set of xy generating vertices*/
    const std::vector<std::pair<TDD_real_t, TDD_real_t> >  xyVertices() const;

    /**This method returns the halflength in local z*/
    double halflengthZ() const;

    /**This method returns the transcript into combined volume*/
    const Acts::Volume* combinedVolume() const;

    /**This method returns the volume envelope*/
    const Acts::Volume* envelope() const;

    /** Output Method for std::ostream */
    std::ostream& dump(std::ostream& sl) const override;

  private:
    template <class T> T& dumpT(T& dt) const;


    void processSubVols() const;
    Acts::PlaneSurface* sideSurf(Transform3D,unsigned int,unsigned int) const;
    bool Xor(bool x, bool y) const;

    bool Left(std::pair<TDD_real_t,TDD_real_t> a, std::pair<TDD_real_t,TDD_real_t> b, std::pair<TDD_real_t,TDD_real_t> c) const;

    bool Intersect(std::pair<TDD_real_t,TDD_real_t> a, std::pair<TDD_real_t,TDD_real_t> b,
		   std::pair<TDD_real_t,TDD_real_t> c, std::pair<TDD_real_t,TDD_real_t> d) const;

    bool InCone(int i, int j, std::vector<std::pair<TDD_real_t,TDD_real_t> > inputVertices) const;

    bool Diagonalie(int  i , int j  ,std::vector<std::pair<TDD_real_t,TDD_real_t> > inputVertices) const;

    bool Diagonal(int i, int j, std::vector<std::pair<TDD_real_t,TDD_real_t> > inputVertices) const;


    std::vector<std::pair<TDD_real_t,TDD_real_t> > TriangulatePolygon(const std::vector<std::pair<TDD_real_t,TDD_real_t> >& Vertices ) const;

    std::vector<std::pair<TDD_real_t,TDD_real_t> > TriangulatePolygonCheck(const std::vector<std::pair<TDD_real_t,TDD_real_t> >& Vertices ) const;

    mutable std::vector<std::pair<double,double> > m_xyVtx; //!< generating xy vertices
    double m_halfX;                                         //!< halflength in x - to define enclosing rectangle
    double m_halfY;                                         //!< halflength in y - to define enclosing rectangle
    double m_halfZ;                                         //!< halflength in z

    mutable int m_ordering;                                 //!< -1 not set/ 1 anticlockwise / 0  clockwise
    mutable const Acts::Volume* m_combinedVolume;            //!< triangulated polygon
    mutable const Acts::Volume* m_envelope;                  //!< simplified envelope

 };

 inline SimplePolygonBrepVolumeBounds* SimplePolygonBrepVolumeBounds::clone() const
 { return new SimplePolygonBrepVolumeBounds(*this); }

 inline const std::vector<std::pair<TDD_real_t,TDD_real_t> > SimplePolygonBrepVolumeBounds::xyVertices() const { return m_xyVtx; }

 inline double SimplePolygonBrepVolumeBounds::halflengthZ() const { return m_halfZ; }

 inline const Acts::Volume* SimplePolygonBrepVolumeBounds::combinedVolume() const { return m_combinedVolume; }

 inline const Acts::Volume* SimplePolygonBrepVolumeBounds::envelope() const { return m_envelope; }

 template <class T> T& SimplePolygonBrepVolumeBounds::dumpT(T& dt) const
 {
     dt << std::setiosflags(std::ios::fixed);
     dt << std::setprecision(7);
     dt << "Acts::SimplePolygonBrepVolumeBounds: (halfZ, xy vertices) = ";
     dt << "( " << m_halfZ << ")";
     for (unsigned int i=0;i<m_xyVtx.size();i++)
       dt << "(" << m_xyVtx[i].first << ","<<m_xyVtx[i].second <<")";
     return dt;
 }


}

#endif // ACTS_VOLUMES_SIMPLEPOLYGONBREPVOLUMEBOUNDS_H

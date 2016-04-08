///////////////////////////////////////////////////////////////////
// CylinderBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACESNOBOUNDS_H
#define ACTS_SURFACESNOBOUNDS_H
    
// Geometry module    
#include "Surfaces/SurfaceBounds.h"
// EventData module
#include "CoreUtils/ParameterDefinitions.h"
// Core module
#include "Algebra/AlgebraDefinitions.h"

class MsgStream;

namespace Acts {
   
  /**
   @class NoBounds
    
   Bounds object for a boundless surface (...)
           
   @author Andreas.Salzburger@cern.ch */
   
   class NoBounds : public SurfaceBounds {
     public:
       /** Default Constructor*/
       NoBounds(){}
       
       /** Destructor */
       ~NoBounds(){}
       
       /** Equality operator */
       virtual bool operator==(const SurfaceBounds& sbo) const override;
       
       /** Return SurfaceBounds for persistency */
       virtual SurfaceBounds::BoundsType type() const final { return SurfaceBounds::Other; }
       
       /** Method inside() returns true for any case */
       virtual bool inside(const Vector2D& locpo, double tol1=0., double tol2=0.) const final;
       virtual bool inside(const Vector2D& locpo, const BoundaryCheck& bchk) const final;
       
       /** This method checks inside bounds in loc1
         - loc1/loc2 correspond to the natural coordinates of the surface */
       virtual bool insideLoc1(const Vector2D& locpo, double tol1=0.) const final;
                        
       /** This method checks inside bounds in loc2 
         - loc1/loc2 correspond to the natural coordinates of the surface */
       virtual bool insideLoc2(const Vector2D& locpo, double tol2=0.) const final;

       /** Minimal distance to boundary (=0 if inside) */
       virtual double minDistance(const Vector2D& pos) const final;
      
       /** Clone method to complete inherited interface */
       virtual NoBounds* clone() const final;

       /** r() method to complete inherited interface */
       virtual double r() const final;
     
       /** Output Method for MsgStream*/
       virtual MsgStream& dump(MsgStream& sl) const final;

       /** Output Method for std::ostream */
       virtual std::ostream& dump(std::ostream& sl) const final;
  
  };

  inline bool NoBounds::operator==(const SurfaceBounds&) const 
    { return true; } 
   
  inline bool NoBounds::inside(const Vector2D&, double, double) const
    { return true; }
	
  inline bool NoBounds::inside(const Vector2D&, const BoundaryCheck&) const
  {  return true;  }

  inline bool NoBounds::insideLoc1(const Vector2D&, double) const
    { return true; }

  inline bool NoBounds::insideLoc2(const Vector2D&, double) const
    { return true; }

  inline double NoBounds::minDistance(const Vector2D&) const
    { return 0.; }

  inline NoBounds* NoBounds::clone() const 
    { return new NoBounds(); }

  inline double NoBounds::r() const 
    { return 0.; }
   
  inline MsgStream& NoBounds::dump(MsgStream& sl) const
  {
       sl << "Acts::NoBounds ... boundless surface" << endreq;
       return sl;
  }
   
  inline std::ostream& NoBounds::dump(std::ostream& sl) const
  {
      sl << "Acts::NoBounds ... boundless surface" << std::endl;
      return sl;   
  }
  
} // end of namespace

#endif // ACTS_SURFACES_NOBOUNDS_H

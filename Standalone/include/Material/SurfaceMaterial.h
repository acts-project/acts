///////////////////////////////////////////////////////////////////
// SurfaceMaterial.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_MATERIAL_SURFACEMATERIAL_H
#define ACTS_MATERIAL_SURFACEMATERIAL_H

// Core module
#include "Core/AlgebraDefinitions.h"
// Geometry module
#include "Material/MaterialProperties.h"
// EventData module
#include "Core/PropDirection.h"
// STD/STL
#include <vector>
#include <memory>

namespace Acts {

  class BinUtility;
  class ElementTable;
    
 /**
   @enum MaterialConcentration

   Simple enum to identify when a material update on a non-structured layer should be done,
   options are alongPre and oppositePre.

  */
  
  enum MaterialConcentration { alongPre  = 1, split =  0, oppositePre = -1 };

  /** 
   @class SurfaceMaterial

   MaterialProperties that are associated to a surface,
   extended by certain special representations
  
   The SurfaceMaterial class inherits from GeometryID,
   in order to allow storing the material in a file and assigning it uniquely.
      
   @author Andreas.Salzburger@cern.ch 
   */

  class SurfaceMaterial {
    
    public:
        
      /**Constructor*/
      SurfaceMaterial() :
       m_splitFactor(1.)
     {}

      /**Constructor*/
      SurfaceMaterial(double splitFactor) :
       m_splitFactor(splitFactor)
      {}

      /**Destructor*/
      virtual ~SurfaceMaterial(){}
      
      /**Pseudo-Constructor clone()*/ 
      virtual SurfaceMaterial* clone() const = 0;
      
      /** Scale operator */
      virtual SurfaceMaterial& operator*=(double scale) = 0;

      /** Return method for full material description of the Surface - from local coordinates */
      virtual const MaterialProperties* material(const Vector2D& lp) const = 0;
      
      /** Return method for full material description of the Surface - from the global coordinates */
      virtual const MaterialProperties* material(const Vector3D& gp) const = 0;
      
      /**Direct access via bins to the MaterialProperties */
      virtual const MaterialProperties* material(size_t ib0, size_t ib1) const = 0;
      
      /** Update the ElementTable */
      void updateElementTable(std::shared_ptr<const ElementTable> ) const { return; }
      
      /** Get the ElementTable */
      const ElementTable* elementTable() const { return nullptr; }
            
      /** Update pre factor */
      double factor(PropDirection pDir, MaterialUpdateStage mStage) const;

      /** Return the BinUtility */
      virtual const BinUtility* binUtility() const = 0;
            
      /** Update the BinUtility if necessary - passing ownership of the utility class*/
      virtual void updateBinning(BinUtility* bu) const = 0;

      /** Output Method for std::ostream, to be overloaded by child classes */
      virtual std::ostream& dump(std::ostream& sl) const = 0;
                                            
    protected :
      double               m_splitFactor;     //!< the split factor in favour of oppositePre

  };

/** inline return methods for the pre/post factors */  
inline double SurfaceMaterial::factor(PropDirection pDir, MaterialUpdateStage mStage ) const 
{ 
    if (mStage == Acts::fullUpdate) return 1.;
    return ( pDir*mStage > 0 ? m_splitFactor : 1.-m_splitFactor ); 
}   

//**Overload of << operator for std::ostream for debug output*/
std::ostream& operator<<( std::ostream& sl, const SurfaceMaterial& sm);
    
} // end of namespace

#endif // ACTS_MATERIAL_SURFACEMATERIAL_H


///////////////////////////////////////////////////////////////////
// HomogeneousSurfaceMaterial.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_MATERIAL_HOMOGENOUSLAYERMATERIAL_H
#define ACTS_MATERIAL_HOMOGENOUSLAYERMATERIAL_H

// Core module
#include "ACTS/Material/SurfaceMaterial.h"
#include "ACTS/Material/MaterialProperties.h"
// STD/STL
#include <vector>
#include "ACTS/Utilities/Definitions.h"

namespace Acts {


  class BinUtility;
  /** 
   @class HomogeneousSurfaceMaterial

   It extends the SurfaceMaterial base class and describes a simple homogeneus material
   descriptions
      
   @author Andreas.Salzburger@cern.ch 
   */

  class HomogeneousSurfaceMaterial : public SurfaceMaterial {
    
    public:
      /** Default Constructor - creates empty HomogeneousSurfaceMaterial */
      HomogeneousSurfaceMaterial();
      
      /**Explizit constructor with only full MaterialProperties, and a split factor*/
      HomogeneousSurfaceMaterial(const MaterialProperties& full, double splitFactor=1.);
      
      /**Copy Constructor */  
      HomogeneousSurfaceMaterial(const HomogeneousSurfaceMaterial& mprop);
      
      /**Destructor*/
      virtual ~HomogeneousSurfaceMaterial();
      
      /**Pseudo-Constructor clone()*/ 
      HomogeneousSurfaceMaterial* clone() const override;
      
      /** Assignment operator */
      HomogeneousSurfaceMaterial& operator=(const HomogeneousSurfaceMaterial& lmp);

      /** Scale operator */
      HomogeneousSurfaceMaterial& operator*=(double scale) override;

      /**Return method for full material description of the Layer - local coordinates*/
      virtual const MaterialProperties* material(const Vector2D& lp) const;
      
      /**Return method for full material description of the Layer - global coordinates*/
      virtual const MaterialProperties* material(const Vector3D& gp) const;

      /**Direct access via bins to the MaterialProperties */
      virtual const MaterialProperties* material(size_t ib0, size_t ib1) const override;
      
      /** Return the BinUtility */
      const BinUtility* binUtility() const  override { return nullptr; }
      
      /** Update the BinUtility if necessary - passing ownership of the utility class*/
      void updateBinning(BinUtility*) const override { }
          
      /** Output Method for std::ostream, to be overloaded by child classes */
      std::ostream& dump(std::ostream& sl) const override;      

    private:
      /** The five different MaterialProperties */
      MaterialProperties*           m_fullMaterial;
                                            
  };
  
inline HomogeneousSurfaceMaterial* HomogeneousSurfaceMaterial::clone() const
  { return new HomogeneousSurfaceMaterial(*this); }  

inline const MaterialProperties* HomogeneousSurfaceMaterial::material(const Vector2D&) const
  { return m_fullMaterial; }
  
inline const MaterialProperties* HomogeneousSurfaceMaterial::material(const Vector3D&) const
  { return m_fullMaterial; }
      
inline const MaterialProperties* HomogeneousSurfaceMaterial::material(size_t, size_t) const
  { return m_fullMaterial; }      
    
} // end of namespace

#endif 



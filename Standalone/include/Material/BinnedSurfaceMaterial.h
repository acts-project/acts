///////////////////////////////////////////////////////////////////
// BinnedSurfaceMaterial.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_MATERIAL_BINNEDSURFACEMATERIAL_H
#define ACTS_MATERIAL_BINNEDSURFACEMATERIAL_H 1

// Geometry module
#include "Material/SurfaceMaterial.h"
#include "Material/MaterialProperties.h"
#include "GeometryUtils/BinUtility.h"
// Core module
#include "Core/AlgebraDefinitions.h"

namespace Acts {


  /** 
   @class BinnedSurfaceMaterial

   It extends the SurfaceMaterial base class and is just an array of it

   @author Andreas.Salzburger@cern.ch 
   */

  class BinnedSurfaceMaterial : public SurfaceMaterial {
    
    public:
      /** Default Constructor - needed by POOL*/    
      BinnedSurfaceMaterial();
  
      /** Default Constructor for emptly material */    
      BinnedSurfaceMaterial(BinUtility& binutility);

      /**Explizit constructor with only full MaterialProperties, 
         and split factors:
          - 1. : oppositePre
          - 0. : alongPre
        ===> 1 Dimensional array

        ATTENTION: Ownership of MaterialProperties objects is given!
       */
      BinnedSurfaceMaterial(const Acts::BinUtility& binutility,
                            const MaterialPropertiesVector& fullProperties,
                            double splitFactor=0.);

      /**Explizit constructor with only full MaterialProperties, 
         and split factors:
          - 1. : oppositePre
          - 0. : alongPre
        ===> 2 Dimensional array

        ATTENTION: Ownership of MaterialProperties objects is given!
       */
      BinnedSurfaceMaterial(const Acts::BinUtility& binutility,
                            const MaterialPropertiesMatrix& fullProperties,
                            double splitFactor=0.);

      /**Copy Constructor */  
      BinnedSurfaceMaterial(const BinnedSurfaceMaterial& mprop);
      
      /**Destructor*/
      virtual ~BinnedSurfaceMaterial();
      
      /**Pseudo-Constructor clone()*/ 
      BinnedSurfaceMaterial* clone() const override;
      
      /** Assignment operator */
      BinnedSurfaceMaterial& operator=(const BinnedSurfaceMaterial& lmp);

      /** Scale operator */
      BinnedSurfaceMaterial& operator*=(double scale) override;

      /** Return the BinUtility */
      const BinUtility* binUtility() const override;
       
      /** Update the BinUtility if necessary - passing ownership of the utility class*/
      void updateBinning(BinUtility* bu) const override; 
       
      /**Return method for full material description of the Layer - for all bins*/
      const MaterialPropertiesMatrix& fullMaterial() const;

      /**Return method for full material description of the Layer - local coordinates */
      const MaterialProperties* material(const Vector2D& lp) const override;
 
      /**Return method for full material description of the Layer - global coordinates */
      const MaterialProperties* material(const Vector3D& gp) const override;
            
      /** Access the single bin */
      const MaterialProperties* material(size_t bin0, size_t bin1 ) const override;
            
      /** Output Method for std::ostream, to be overloaded by child classes */
      std::ostream& dump(std::ostream& sl) const override;      

    private:

      mutable BinUtility*       m_binUtility; //!< the helper for the bin finding
 
      /** The five different MaterialProperties */
      MaterialPropertiesMatrix m_fullMaterial;

      /** helper method - to clear the material*/
      void clearMaterial();
                                            
      /** helper method - to refill the material  */
      void fillMaterial(const MaterialPropertiesMatrix& matMatrix);

  };


inline const BinUtility* BinnedSurfaceMaterial::binUtility() const
  { return m_binUtility; }
  
  inline const MaterialPropertiesMatrix& BinnedSurfaceMaterial::fullMaterial() const
  { return m_fullMaterial; }
  
  inline const MaterialProperties* BinnedSurfaceMaterial::material(size_t bin0, size_t bin1 ) const 
  {
     return m_fullMaterial[bin1][bin0];
  }
  
  inline void BinnedSurfaceMaterial::updateBinning(BinUtility* bu) const {
      if (bu){
          delete m_binUtility;
          m_binUtility = bu;
      }
  }


}

#endif

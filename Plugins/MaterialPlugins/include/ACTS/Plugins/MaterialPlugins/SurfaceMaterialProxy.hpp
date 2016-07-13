// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceMaterialProxy.hpp, ACTS project MaterialPlugins
///////////////////////////////////////////////////////////////////

#ifndef ACTS_MATERIAL_SURFACEMATERIALPROXY_H
#define ACTS_MATERIAL_SURFACEMATERIALPROXY_H

#include "ACTS/Material/SurfaceMaterial.hpp"
#include "ACTS/Utilities/BinUtility.hpp"

namespace Acts {
    
    /// @class SurfaceMaterialProxy
    ///
    /// @brief proxy to SurfaceMaterial hand over BinUtility
    ///
    /// The SurfaceMaterialProxy class acts as a proxy to the SurfaceMaterial
    /// to mark the layers and surfaces on which the material should be mapped on
    /// at construction time of the geometry and to hand over the granularity of
    /// of the material map with the bin utility.
    
    class SurfaceMaterialProxy : public SurfaceMaterial {
        
    public:
        /// Constructor with BinUtility
        /// @param binutility two dimensional BinUtility determining the granularity and binning of the material on the surface/layer
        SurfaceMaterialProxy(BinUtility& binutility);
        /// Copy constuctor
        SurfaceMaterialProxy(const SurfaceMaterialProxy& smproxy);
        /// Destructor
        virtual ~SurfaceMaterialProxy() = default;
        /// Implicit contructor
        /// - uses the copy constructor
        SurfaceMaterialProxy* clone() const override;
        /// Scale operator
        virtual SurfaceMaterial& operator*=(double scale) override;
        
        /// Return method for full material description of the Surface - from local coordinates
        virtual const MaterialProperties* material(const Vector2D& lp) const override;
        
        /// Return method for full material description of the Surface - from the global coordinates
        virtual const MaterialProperties* material(const Vector3D& gp) const override;
        
        /// Direct access via bins to the MaterialProperties
        virtual const MaterialProperties* material(size_t ib0, size_t ib1) const override;
        
        
        /// Return the BinUtility
        virtual const BinUtility* binUtility() const override;
        
        /// Update the BinUtility if necessary - passing ownership of the utility class
        virtual void updateBinning(BinUtility* bu) override;
        
        /// Output Method for std::ostream, to be overloaded by child classes
        virtual std::ostream& dump(std::ostream& sl) const override;
        
    private:
        /// two dimensional BinUtility determining the granularity and binning of the material on the surface/layer
        BinUtility* m_binUtility;
        
    };
    
}

inline const Acts::MaterialProperties* Acts::SurfaceMaterialProxy::material(const Vector2D&) const {return nullptr;}

inline const Acts::MaterialProperties* Acts::SurfaceMaterialProxy::material(const Vector3D&) const {return nullptr;}

inline const Acts::MaterialProperties* Acts::SurfaceMaterialProxy::material(size_t, size_t) const {return nullptr;}

inline const Acts::BinUtility* Acts::SurfaceMaterialProxy::binUtility() const {return m_binUtility;}

inline void Acts::SurfaceMaterialProxy::updateBinning(BinUtility* bu)
{
    if (bu) {
        delete m_binUtility;
        m_binUtility = bu;
    }
}



#endif // ACTS_MATERIAL_SURFACEMATERIALPROXY_H
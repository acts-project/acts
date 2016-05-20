// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DetectorElementBase.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DETELEMENTBASE_DETELEMENTBASE_H
#define ACTS_DETELEMENTBASE_DETELEMENTBASE_H 1

#ifdef ACTS_GEOMETRY_DETELEMENT_PLUGIN
#include ACTS_GEOMETRY_DETELEMENT_PLUGIN
#else

// Algebra and Identifier
#include "ACTS/Utilities/Identifier.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include <memory>
#include <vector>

namespace Acts
{

    class Surface;
    class SurfaceBounds;

    /** @class DetectorElementBase

     This is the base class for all tracking detector elements
     with read-out relevant information.
     
     If a DetectorElement has a second element (or even a triple setup) 
     that would naturally fall into the same bin, one can register that as a binmember.
    
     DetectorElements close by can be registered as neighbours as this will help
     the navigation.
     
     @author Andreas.Salzburger@cern.ch

     */

    class DetectorElementBase {

      public:

        /** Constructor */
        DetectorElementBase(){}

        /** virtual Destructor */
        virtual ~DetectorElementBase(){}

        /** Identifier */
        virtual Identifier identify() const = 0;
                
        /** Return local to global transform (optionally associated with an identifier) */
        virtual const Transform3D& transform(const Identifier& identifier = Identifier()) const = 0;
                
        /** Return the surface associated with this detector element (optionally associated with an identifier) */
        virtual const Surface& surface(const Identifier& identifier = Identifier()) const = 0;
        
        /** Return the surface bounds (optionally with an associated detector element) */
        virtual const SurfaceBounds& bounds(const Identifier& identifier = Identifier()) const = 0;

        /** Returns the full list of all detection surfaces associated to this detector element */
        virtual const std::vector< std::shared_ptr<const Surface> >& surfaces() const = 0;
                
        /** Return the center of the element (optionally associated with an identifier) */
        virtual const Vector3D& center(const Identifier& identifier = Identifier()) const = 0;
        
        /** Return the normal of the element (optionally associated with an identifier) */
        virtual const Vector3D& normal(const Identifier& identifier = Identifier()) const = 0;
        
        /** Returns the thickness of the module */
        virtual double thickness() const = 0;

        /** Bin members for fast access */
        const std::vector<const DetectorElementBase*>& binmembers() const;
        
        /** Neighbours for fast access */
        void registerBinmembers(std::vector<const DetectorElementBase*>& binmembers) const;
                
        /** Neighbours for fast access */
        const std::vector<const DetectorElementBase*>& neighbours() const;
        
        /** Neighbours for fast access */
        void registerNeighbours(std::vector<const DetectorElementBase*>& neighbours) const;
        
    private:
        mutable std::vector<const DetectorElementBase*>    m_binmembers;
        mutable std::vector<const DetectorElementBase*>    m_neighbours;
        
    };
    
    inline const std::vector<const DetectorElementBase*>& DetectorElementBase::binmembers() const { return m_binmembers; } 
       
    inline const std::vector<const DetectorElementBase*>& DetectorElementBase::neighbours() const { return m_neighbours; }    
    

    inline void DetectorElementBase::registerBinmembers(std::vector<const DetectorElementBase*>& binmembers) const
    { 
        for (auto& bmember : binmembers){
          // only fill if it's not yet registered    
          if (find(m_binmembers.begin(),m_binmembers.end(),bmember) == m_binmembers.end())
               m_binmembers.push_back(bmember);
      }
    }

    inline void DetectorElementBase::registerNeighbours(std::vector<const DetectorElementBase*>& neighbours) const
    { 
        for (auto& neighbour : neighbours){
          // only fill if it's not yet registered    
          if (find(m_neighbours.begin(),m_neighbours.end(),neighbour) == m_neighbours.end())
               m_neighbours.push_back(neighbour);
      }
    }
    
    
}//end of ns
#endif

#endif // ACTS_GEOMETRY_DETELEMENT_PLUGIN

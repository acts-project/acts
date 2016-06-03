// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BinUtility.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYUTILS_BINUTILITY_H
#define ACTS_GEOMETRYUTILS_BINUTILITY_H 1

// Core module
#include "ACTS/Utilities/BinningData.hpp"
#include "ACTS/Utilities/BinningType.hpp"
// STL
#include <vector>
#include <memory>
#include "Definitions.hpp"

//debug
#include <iostream>

namespace Acts {

/** @class BinUtility
    
    The BinUtility class that translated global and local position into a bins of a BinnedArray,
    most performant is equidistant binning without a transform, however, 
    optionally a transform can be provided, e.g. for binning on shifted object,
    the transform is usually shared with the geometric object the Array is defined on,
    for performance reasons, also the inverse transform is stored.
        
  */

   class BinUtility {

      public:
       /** Constructor for equidistant */
        BinUtility() :
         m_binningData(),
         m_transform(nullptr),
         m_itransform(nullptr)
         {
           m_binningData.reserve(3);
         }

       /** Constructor from BinningData directly */
        BinUtility(const BinningData& bData, std::shared_ptr<Transform3D> tForm=nullptr) :
         m_binningData(),
         m_transform(tForm),
         m_itransform(tForm ? new Transform3D(tForm->inverse()) : nullptr)
         {
           m_binningData.reserve(3);
           m_binningData.push_back(bData);
         }

        /** Constructor for equidistant  */
        BinUtility(size_t bins, float min, float max, BinningOption opt = open, BinningValue value = binX, std::shared_ptr<Transform3D> tForm=nullptr) :
          m_binningData(),
          m_transform(tForm),
          m_itransform(tForm ? new Transform3D(tForm->inverse()) : nullptr)
          {
            m_binningData.reserve(3);
            m_binningData.push_back(BinningData(opt, value, bins, min, max));
          }

        /** Constructor for arbitrary */
        BinUtility(std::vector<float>& bValues, BinningOption opt = open, BinningValue value = binPhi, std::shared_ptr<Transform3D> tForm=nullptr) :
	      m_binningData(),
          m_transform(tForm),
         m_itransform(tForm ? new Transform3D(tForm->inverse()) : nullptr)
        {
	        m_binningData.reserve(3);
	        m_binningData.push_back(BinningData(opt, value, bValues));
	    }

	    /** Constructor for binH */
	    BinUtility(float phiRef, std::vector<std::pair<int,float> >& bValues, std::shared_ptr<Transform3D> tForm=nullptr) :
	      m_binningData(),
          m_transform(tForm),
          m_itransform(tForm ? new Transform3D(tForm->inverse()) : nullptr)
        {
            m_binningData.reserve(3);
            m_binningData.push_back(BinningData(open,phiRef,bValues));
        }

        /** Copy constructor */
        BinUtility(const BinUtility& sbu ) :
         m_binningData(sbu.m_binningData),
         m_transform(sbu.m_transform),
         m_itransform(sbu.m_transform ? new Transform3D(sbu.m_transform->inverse()) : nullptr)
       {}

        /** Assignment operator Constructor */
        BinUtility& operator=(const BinUtility& sbu ){
            if ( this != &sbu ){
                m_binningData = sbu.m_binningData;
                m_transform   = sbu.m_transform;
                m_itransform  = sbu.m_transform ? std::unique_ptr<Transform3D>(new Transform3D(sbu.m_transform->inverse())) : nullptr;
            }
            return (*this);
        }

        /** Operator++ to make multidimensional BinUtility */
        BinUtility& operator+= ( const BinUtility& gbu) throw (std::string) {
            const std::vector<BinningData>& bData = gbu.binningData();
            if (m_binningData.size() + bData.size() > 3)
                throw "BinUtility does not support dim > 3";
            m_binningData.insert(m_binningData.end(), bData.begin(), bData.end());
            return (*this);
        }

        /** Virtual Destructor */
        virtual ~BinUtility(){}

        /** Implizit Constructor */
        virtual BinUtility* clone() const { return new BinUtility(*this); }

        /** return the binning data */
        const std::vector<BinningData>& binningData() const { return m_binningData; }
          
        /** Bin from a 3D vector (already in binning frame) - optionally the itransform is applied */
        size_t bin(const Vector3D& position, size_t ba=0) const throw (std::string)
        {
          if (ba >= m_binningData.size()) 
                throw "dimension out of bounds"; 
          size_t bEval = m_itransform ? m_binningData[ba].searchGlobal((*m_itransform)*position) : m_binningData[ba].searchGlobal(position);          
          return ( bEval > bins(ba)-1 ? bins(ba)-1 : bEval );     //!< @TODO ST additional protection : DEBUG source
        }
        
        /** Bin neighbour range */
        std::vector<size_t> neighbourRange(const Vector3D& position, size_t ba=0) const 
        {
            std::vector<size_t> neighbourRange;
            size_t cbin  = bin(position,ba);
            size_t pbin = cbin ? cbin-1 : ( ( m_binningData[ba].option == open ) ? cbin : m_binningData[ba].bins()-1 );
            size_t nbin = (cbin < m_binningData[ba].bins()-1) ? cbin+1 : 0;
            if (pbin != cbin ) neighbourRange.push_back(pbin);
            neighbourRange.push_back(cbin);
            if (nbin != cbin ) neighbourRange.push_back(nbin);
            return neighbourRange;
        }    
        
        /** Bin from a 3D vector (already in binning frame) */
        size_t entry(const Vector3D& position, size_t ba=0) const throw (std::string)
        { if (ba >= m_binningData.size()) 
	    throw "dimension out of bounds"; 
          return (m_itransform ?  m_binningData[ba].entry((*m_itransform)*position) : m_binningData[ba].entry(position));
        }       
        
        /** Bin from a 3D vector (already in binning frame) */
        size_t next(const Vector3D& position, const Vector3D& direction, size_t ba=0) const throw (std::string)
        { if (ba >= m_binningData.size()) 
                throw "dimension out of bounds"; 
          return (m_itransform ? m_binningData[ba].next((*m_itransform)*position, (m_itransform->linear())*direction) : m_binningData[ba].next(position, direction));
        }       

        /** Return the oder direciton for fast interlinking */
        int nextDirection(const Vector3D& position, const Vector3D& direction, size_t ba=0) const throw (std::string)  {
             if (ba >= m_binningData.size())
                    throw "dimension out of bounds";
              return m_binningData[ba].nextDirection(position, direction);
        }

        /** Distance estimate to next bin  */
        std::pair<size_t,float> distanceToNext(const Vector3D& position, const Vector3D& direction, size_t ba=0) const throw (std::string)
        { if (ba >= m_binningData.size()) 
                throw "dimension out of bounds"; 
          return (m_itransform ?  m_binningData[ba].distanceToNext((*m_itransform)*position, (m_itransform->linear())*direction) : m_binningData[ba].distanceToNext(position, direction));
        }       
                
        /** Bin from a 2D vector (following local parameters defintitions) - no optional transform applied 
            - USE WITH CARE !!
              You need to check if your local position is actually in the binning frame of the BinUtility */
        size_t bin(const Vector2D& lposition, size_t ba=0) const throw (std::string) 
        {    
            if (ba >= m_binningData.size()) 
                    throw "dimension out of bounds"; 
              return m_binningData[ba].searchLocal(lposition);
        }
        
        /** Check if bin is inside from Vector3D - optional transform applied */
        bool inside(const Vector3D& position ) const {
	        std::vector<BinningData>::const_iterator bdIter = m_binningData.begin();
	        for ( ; bdIter != m_binningData.end(); ++bdIter) {
              if (m_itransform && !(*bdIter).inside((*m_itransform)*position)) return false;    
	          else if ( !(*bdIter).inside(position)) return false;
		}
	        return true;
        }
        
        /** Check if bin is inside from Vector2D - no optional transform applied */
        bool inside(const Vector2D& lposition) const {
            return true;
            std::vector<BinningData>::const_iterator bdIter =  m_binningData.begin();
            for ( ; bdIter != m_binningData.end(); ++bdIter)
              if ( !(*bdIter).inside(lposition)) return false;
            return true;
        }

        /** First bin maximal value */
        size_t dimensions() const {
	        return m_binningData.size();
        }

        /** First bin maximal value */
        size_t max(size_t ba=0) const {
	        if (ba >= m_binningData.size()) return 0;
	        return (m_binningData[ba].bins()-1);
        }

        /** Number of bins */
        size_t bins(size_t ba=0) const {
	        if (ba >= m_binningData.size()) return 0;
	        return (m_binningData[ba].bins());
        }

        /** The type/value of the binning */
        BinningValue binningValue(size_t ba=0) const throw (std::string) {
            if (ba >= m_binningData.size())
                throw "dimension out of bounds";
            return (m_binningData[ba].binvalue);
         }

        /** bin->BinningValue navigation : pos=+-1. edges/ 0. bin center */
        float binPosition( size_t bin, float pos, size_t ba=0 ) const {
            if (ba >= m_binningData.size())
                throw "dimension out of bounds";
            return (m_binningData[ba].binPosition( bin, pos ));
	}

        /** Clear the data. */
        void clear() { m_binningData.clear(); }

       /** Output Method for std::ostream, to be overloaded by child classes */
       std::ostream& dump(std::ostream& sl) const {
            sl << "BinUtility for " << m_binningData.size() << "- dimensional array:" << std::endl;
            std::vector<BinningData>::const_iterator bdIter = m_binningData.begin();
            for (size_t ibd = 0 ; bdIter != m_binningData.end(); ++bdIter, ++ibd ){
                sl << "dimension     : " << ibd                                          << std::endl;
                sl << " - type       : " << size_t((*bdIter).type)                       << std::endl;
                sl << " - option     : " << size_t((*bdIter).option)                     << std::endl;
                sl << " - value      : " << size_t((*bdIter).binvalue)                   << std::endl;
                sl << " - bins       : " << (*bdIter).bins()                             << std::endl;
                sl << " - min/max    : " << (*bdIter).min << " / " << (*bdIter).max      << std::endl;
                if ((*bdIter).type == equidistant)
                sl << " - step       : " << (*bdIter).step                               << std::endl;
                sl << " - boundaries : | ";
                std::vector<float>::const_iterator bIter = (*bdIter).boundaries().begin();
                for ( ; bIter != (*bdIter).boundaries().end(); ++bIter )
                    sl << (*bIter) << " | ";
                sl << std::endl;
            }
            return sl;
       }

     private :
        std::vector<BinningData>      m_binningData;
        std::shared_ptr<Transform3D>  m_transform;
        std::unique_ptr<Transform3D>  m_itransform;          
            
  };    

/**Overload of << operator for std::ostream for debug output*/
std::ostream& operator << ( std::ostream& sl, const BinUtility& bgen);

} // end of namespace Acts

#endif // ACTS_GEOMETRYUTILS_BINUTILITY_H


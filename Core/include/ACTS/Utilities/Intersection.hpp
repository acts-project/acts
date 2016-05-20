// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// Intersection.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYUTILS_INTERSECTION_H
#define ACTS_GEOMETRYUTILS_INTERSECTION_H 1

#include "Definitions.hpp"

namespace Acts {
       
   /**
     @struct Intersection 
     
     */
      
    struct Intersection {
        
       Vector3D      position;
       double        pathLength;
       double        distance;
       bool          valid;
        
       Intersection(const Vector3D& sinter,
                     double slenght,
                     bool svalid,
                     double dist = 0.) :
         position(sinter),
         pathLength(slenght),
         distance(dist),
         valid(svalid)
       {}
          
       Intersection() :
         position(Vector3D(0.,0.,0.)),
         pathLength(0.),
         distance(0.),
         valid(false)
       {}  
          
       // smaller operator for sorting
       bool operator< (const Intersection& si ) const 
           { return (valid && pathLength < si.pathLength); }
    };


   /** class extensions to return also the object */
   template <class T> class ObjectIntersection {
     public:  
       Intersection  intersection;
       mutable const T*    object;
       int                 pDirection;
   
       /** Default constructor */
       ObjectIntersection():
         intersection(),
         object(nullptr),
         pDirection(0)
       {}
           
       /** Object intersection */
       ObjectIntersection(const Intersection& sInter,
                          const T* sObject,
                          int dir = 1):
         intersection(sInter),
         object(sObject),
         pDirection(dir)
       {}
         
        /** smaller operator for ordering & sorting */
        bool operator< ( const ObjectIntersection<T>& oi ) const {
             return ( intersection <  oi.intersection );
        }
    };

    /** Class extension to return the object, a represenation & the result */
    template <class T, class R, class S> class FullIntersection {
      public:  
        Intersection        intersection;
        mutable const T*    object;
        mutable const R*    representation;
        mutable const S*    result;
        int                 pDirection;
    
        /** Full intersection */
        FullIntersection(const Intersection& sInter,
                         const T* sObject,
                         const R* sRepresentation,
                         const S* sResult,
                         int   dir=1):
          intersection(sInter),
          object(sObject),
          representation(sRepresentation),
          result(sResult),
          pDirection(dir)
        {}
          
         /** smaller operator for ordering & sorting */
         bool operator< ( const FullIntersection<T,R, S>& oi ) const {
              return ( intersection <  oi.intersection );
         }
    };

}

#endif // ACTS_GEOMETRYUTILS_INTERSECTION_H 

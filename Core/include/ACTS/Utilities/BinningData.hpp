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
#ifndef ACTS_GEOMETRYUTILS_BINNINGDATA_H
#define ACTS_GEOMETRYUTILS_BINNINGDATA_H 1

// Core module
#include "ACTS/Utilities/BinningType.hpp"

#include <vector>
#include <utility>
#include <cmath>
#include "Definitions.hpp"
//debug
#include <iostream>

namespace Acts {

   /** @class BinningData

        This class holds all the data necessary for the bin calculation

        phi has a very particular behaviour:
        - there's the change around +/- PI

        @TODO add description of convention for sub structure - it can be multiplicative or additive
        multiplicative : each bin has the same sub structure, i.e. first binnning structure is equidistant
        additive : sub structure replaces one bin (and one bin only)

   */

   class BinningData {

     public :
       /** holding all the data for binning calculatuion */
       BinningType                         type;
       BinningOption                       option;
       BinningValue                        binvalue;
       float                               min;
       float                               max;
       float                               step;
       // specialized for H binning @TODO write documetnation
       float                               refphi;    // ref Phi for binH
       std::vector<std::pair<int,float> >  hbounds;   // boundary for binH
       // sub structure
       BinningData*                        subBinningData;      // describe some sub binning structure
       bool                                subBinningAdditive;  // sub binnint is either additive or multipicative

       /** Constructor for equidistant binning - and optional sub structure can be mulitplicative or additive */
       BinningData ( BinningOption bOption,
                     BinningValue bValue,
                     size_t bBins,
                     float  bMin,
                     float  bMax,
                     BinningData* sBinData = nullptr,
                     bool sBinAdditive = false) :
	    type(equidistant),
	    option(bOption),
	    binvalue(bValue),
	    min(bMin),
	    max(bMax),
	    step( (bMax-bMin)/bBins ),
  	    refphi(0.),
  	    hbounds(std::vector<std::pair<int,float> >()),
        subBinningData( sBinData ),
        subBinningAdditive(sBinAdditive),
	    m_bins(bBins),
  	    m_boundaries(std::vector<float>()),
	    m_totalBins(bBins),
  	    m_totalBoundaries(std::vector<float>()),
        m_functionPtr(nullptr),
        m_mixPtr(nullptr)
	 {
        // set to equidistant search
        m_functionPtr = &searchEquidstantWithBoundary;
        // fill the boundary vector for fast access to center & boundaries
        m_boundaries.reserve(m_bins+1);
        for (size_t ib=0; ib < m_bins+1; ++ib)
             m_boundaries.push_back(min+ib*step);
        // the binning data has sub structure - multiplicative or additive
        checkSubStructure();
     }

     /** Constructor for equidistant binning - and optional sub structure can only be additive */
     BinningData ( BinningOption bOption,
                   BinningValue bValue,
                   const std::vector<float> bBoundaries,
                   BinningData* sBinData = nullptr) :
        type(arbitrary),
        option(bOption),
        binvalue(bValue),
        min(0.),
        max(0.),
        step(0.),
	    refphi(0.),
	    hbounds(std::vector<std::pair<int,float> >()),
        subBinningData( sBinData ),
        subBinningAdditive(true),
        m_bins(bBoundaries.size()-1),
	    m_boundaries(bBoundaries),
        m_totalBins(bBoundaries.size()),
	    m_totalBoundaries(bBoundaries),
        m_functionPtr(nullptr),
        m_mixPtr(nullptr)
    {
        // assert a no-size case
        assert(m_boundaries.size()>1);
        min = m_boundaries[0];
        max = m_boundaries[m_boundaries.size()-1];
        // set to equidistant search
        m_functionPtr = m_bins < 50 ? &searchInVectorWithBoundary : &binarySearchWithBoundary;
        // the binning data has sub structure - multiplicative
        checkSubStructure();
     }

     /** Constructor for binH type : non-equidistant binning assumed */
     BinningData ( BinningOption bOption, float  bRefPhi,
                   const std::vector< std::pair<int,float> >& bBoundaries) :
       type(arbitrary),
       option(bOption),
       binvalue(binH),
       min(bBoundaries.front().second),
       max(bBoundaries.back().second),
       step(1.),                        // non-zero value needed for next()
       refphi(bRefPhi),
       hbounds(bBoundaries),
       subBinningData(nullptr),
       subBinningAdditive(false),
       m_bins(bOption==open? bBoundaries.size()-1 : bBoundaries.size()),
       m_boundaries(std::vector<float>()),
       m_totalBins(0),
 	   m_totalBoundaries(std::vector<float>()),
       m_functionPtr( nullptr ),
       m_mixPtr( &searchInVectorWithMixedBoundary )
    {}

    /** Copy constructor */
    BinningData (const BinningData& bdata) :
      type(bdata.type),
      option(bdata.option),
      binvalue(bdata.binvalue),
      min(bdata.min),
      max(bdata.max),
      step(bdata.step),                        // non-zero value needed for next()
      refphi(bdata.refphi),
      hbounds(bdata.hbounds),
      subBinningData ( nullptr ),
      subBinningAdditive(bdata.subBinningAdditive),
      m_bins(bdata.m_bins),
      m_boundaries(bdata.m_boundaries),
      m_totalBins(bdata.m_totalBins),
      m_totalBoundaries(bdata.m_totalBoundaries),
      m_functionPtr( nullptr ),
      m_mixPtr( nullptr )
     {
        // get the binning data
        subBinningData = bdata.subBinningData ? new BinningData(*bdata.subBinningData) : nullptr;
        // set the pointer depending on the type
        if (binvalue == binH){
            // set the mixed function ptr
            m_mixPtr =  &searchInVectorWithMixedBoundary;
        } else {
            // set the correct function pointer
            if (type == equidistant) m_functionPtr = &searchEquidstantWithBoundary;
            else m_functionPtr = m_bins < 50 ? &searchInVectorWithBoundary : &binarySearchWithBoundary;
        }
     }

     /** assignment operator */
     BinningData& operator=(const BinningData& bdata)
     {
         if (this!=&bdata){
             type                  = bdata.type;
             option                = bdata.option;
             binvalue              = bdata.binvalue;
             min                   = bdata.min;
             max                   = bdata.max;
             step                  = bdata.step;
             refphi                = bdata.refphi;
             hbounds               = bdata.hbounds;
             subBinningAdditive    = bdata.subBinningAdditive;
             subBinningData        = bdata.subBinningData ? new BinningData(*bdata.subBinningData) : nullptr;
             m_bins                = bdata.m_bins;
             m_boundaries          = bdata.m_boundaries;
             m_totalBins           = bdata.m_totalBins;
             m_totalBoundaries     = bdata.m_totalBoundaries;
             // set the pointer depending on the type
             if (binvalue == binH){
                 // set the mixed function ptr
                 m_mixPtr =  &searchInVectorWithMixedBoundary;
             } else {
                 // set the correct function pointer
                 if (type == equidistant) m_functionPtr = &searchEquidstantWithBoundary;
                 else m_functionPtr = m_bins < 50 ? &searchInVectorWithBoundary : &binarySearchWithBoundary;
             }
         }
         return (*this);
      }

     /** Destructor */
     ~BinningData(){
         delete subBinningData;
     }

     /** return the number of bins - including sub bins */
     size_t bins() const { return m_totalBins; }

     /** return the boundaries  - including sub boundaries */
     const std::vector<float>& boundaries() const {
         if (subBinningData) return m_totalBoundaries;
         return m_boundaries;
     }

     /** take the right float value - assumes the correct local position expression */
     float value(const Vector2D& lposition) const {
       // ordered after occurence
       if ( binvalue == binR || binvalue == binRPhi || binvalue == binX || binvalue == binH ) return lposition[0];
       if ( binvalue == binPhi ) return gaugePhi(lposition[1]);
       return lposition[1];
     }

     /** take the right float value */
     float value(const Vector3D& position) const {
       // ordered after occurence
       if ( binvalue == binR || binvalue == binH ) return (position.perp());
       if ( binvalue == binRPhi ) return  (position.perp()*position.phi());
       if ( binvalue == binEta ) return (position.eta());
       if ( binvalue < 3  ) return (position[binvalue]);
       // phi gauging
       return gaugePhi(position.phi());
     }

     /** gauge phi */
     float gaugePhi(float phi) const {
         if (max > M_PI && phi < min && phi < 0.){
             phi += 2.*M_PI;
         }
         return phi;
     }

     /** take float values for binH */
     std::pair<float,float> valueH(const Vector2D& lposition ) const {
       return (std::pair<double,double> (lposition[0],lposition[0]*cos(fabs(refphi-lposition[1]))));
     }

     /** take float values for binH */
     std::pair<float,float> valueH(const Vector3D& position ) const {
       return (std::pair<double,double> (position.perp(),position.perp()*cos(fabs(position.phi()-refphi))));
     }

     /** Check if bin is inside from Vector3D */
     bool inside(const Vector3D& position) const {
       // closed one is always inside
       if (option == closed) return true;
       // all other options except value H
       if ( binvalue != binH ) {
	    float val = value(position);
	    return ( val > min-0.001 && val < max+0.001);
       }
       // value H case
       std::pair<double,double> valH = valueH(position);
       float valmin = hbounds.front().first==0 ? valH.first : valH.second;
       float valmax = hbounds.back().first==0 ?  valH.first : valH.second;
       return ( valmin > min-0.001 && valmax < max+0.001);
     }

     /** Check if bin is inside from Vector2D */
     bool inside(const Vector2D& lp) const  {
       if (option == closed) return true;
       if ( binvalue != binH ) {
	    float val = value(lp);
	    return ( val > min-0.001 && val < max+0.001);
       }
       std::pair<double,double> valH = valueH(lp);
       float valmin = hbounds.front().first==0 ? valH.first : valH.second;
       float valmax = hbounds.back().first==0 ?  valH.first : valH.second;
       return ( valmin > min-0.001 && valmax < max+0.001);
     }

     /** generic search from a 2D position --- corresponds to local coordinate schema */
     size_t searchLocal(const Vector2D& lposition) const {
       return (binvalue==binH) ? searchH(valueH(lposition)) : search(value(lposition));
     }

     /** generic search from a 3D position */
     size_t searchGlobal(const Vector3D& position) const {
       return (binvalue==binH) ? searchH(valueH(position)) : search(value(position));
     }

     /** generic search - forwards to correct function pointer */
     size_t search(float value) const {
         assert(m_functionPtr != nullptr);
         return (!subBinningData) ? (*m_functionPtr)(value, *this) : searchWithSubStructure(value);
     }

     /** generic search  sub structure - forwards to correct function pointer */
     size_t searchWithSubStructure(float value) const {
         // find the masterbin with the correct function pointer
          size_t masterbin = (*m_functionPtr)(value, *this);
          // additive sub binning -
          if (subBinningAdditive) {
              // no gauging done, for additive sub structure
              return masterbin + subBinningData->search(value);
          }
          // gauge the value to the subBinData
          float gvalue = value - masterbin*(subBinningData->max-subBinningData->min);
          // now go / additive or multiplicative
          size_t subbin = subBinningData->search(gvalue);
          // now return
          return masterbin*subBinningData->bins() + subbin;
     }

     /** generic search - forwards to correct function pointer */
     size_t searchH(std::pair<double,double> value) const { assert(m_mixPtr != nullptr); return (*m_mixPtr)(value, *this); }

     /** the entry bin */
     size_t entry(const Vector3D& position ) const {
       size_t bin = (binvalue==binH) ? searchH(valueH(position)) : search(value(position));
       return (bin < m_bins-bin) ? bin : m_bins-1;
     }

     /** layer next direction is needed  */
     int nextDirection(const Vector3D& position, const Vector3D& dir) const {
           float     val  = (binvalue==binH) ? valueH(position).first : value(position);
           Vector3D probe = position+dir.normalized();
           float  nextval = (binvalue==binH) ? valueH(probe).first : value(probe);
           return (nextval > val ) ? 1 : -1;
     }

     /** the next bin : gives -1 if the next one is outside */
     size_t next(const Vector3D& position, const Vector3D& dir) const {
       float   val    = value(position);
       Vector3D probe = position+0.5*step*dir.normalized();
       float  nextval = value( probe);
       int bin        = (binvalue==binH) ? searchH(valueH(position)) : search(val);
       bin = (nextval > val && bin != int(m_bins-1)) ? bin+1 : (bin) ? bin-1 : 0;
       // closed setup
       if ( option == closed)
           return ( bin < 0 || bin+1 > int(m_bins) ) ? ( (bin < 0 ) ? m_bins-1 : 0  ) : bin;
       // open setup
       return bin;
     }

     /** distance to the next bin : gives -1 if the next one is outside */
     std::pair<size_t, float>  distanceToNext(const Vector3D& position, const Vector3D& dir) const {
         // current value
         float   val    = (binvalue==binH) ? valueH(position).first : value(position);
         // probe value
         Vector3D probe = position+0.5*step*dir.normalized();
         float  nextval = (binvalue==binH) ? valueH(probe).first : value(probe);
         // current bin
         int bin0        = (binvalue==binH) ? searchH(valueH(position)) : search(val);
         // next bin
         int bin =  (nextval-val) >0. ? bin0+1 : bin0-1;
         if (bin > int(m_bins)-1) bin = (option==closed) ? 0 : bin0;
         if (bin < 0) bin = (option==closed) ? m_bins-1 : 0;

         // boundary value
         float bval = 0.;
         if  (binvalue==binH) {
             bval = (nextval >val) ? hbounds[bin0+1].second : hbounds[bin0].second;    // non-equidistant

             // may need to recalculate current value and probe
             if (nextval>val) {
                 if ( hbounds[bin0+1].first>0 ) {
                     val = valueH(position).second;
                     nextval = valueH(probe).second;
                 }
             } else {
                 if ( hbounds[bin0].first>0 ) {
                     val = valueH(position).second;
                     nextval = valueH(probe).second;
                 }
             }
         } else {
             bval = (nextval >val) ? m_boundaries[bin0+1] : m_boundaries[bin0];    // non-equidistant
             if (type == equidistant ) bval = min+bin0*step;
         }
         // differential
         float dd = 2*(nextval-val)/step;
         // distance estimate
         float dist = fabs(dd)>1.e-06 ? (bval-val)/dd : 1.e06 ;
         return std::pair<size_t, float> (bin,dist);
     }

     /** bin->BinningValue navigation : pos=+-1. edges/ 0. bin center */
     float binPosition( size_t bin, float pos ) const {

       float bmin = (binvalue==binH) ? hbounds[bin].second : m_boundaries[bin] ;
       float bmax = (binvalue==binH) ? hbounds[bin+1].second :
             bin+1<m_boundaries.size() ? m_boundaries[bin+1] : m_boundaries[bin]+step  ;

       return ( bmin + 0.5*(pos+1.)*(bmax-bmin) );

     }


    private:
      /** Private members */
      size_t                              m_bins;
      std::vector<float>                  m_boundaries;
      size_t                              m_totalBins;       //!< including potential substructure
      std::vector<float>                  m_totalBoundaries; //!< including potential substructure
      /** the pointer to the function to be used */
      size_t (*m_functionPtr) (float, const BinningData&);
      size_t (*m_mixPtr) (std::pair<float,float>, const BinningData&);

      /** helper method to set the sub structure */
      void checkSubStructure(){
         // sub structure is only checked when sBinData is ndefined
         if (subBinningData){
             m_totalBoundaries.clear();
             // (A) additive sub structure
             if (subBinningAdditive){
                 // one bin is replaced by the sub bins
                 m_totalBins = m_bins+subBinningData->bins()-1;
                 // the tricky one - exchange one bin by many others
                 m_totalBoundaries.reserve(m_totalBins+1);
                 // get the sub bin boundaries
                 const std::vector<float>& subBinBoundaries = subBinningData->boundaries();
                 float sBinMin = subBinBoundaries[0];
                 // get the min value of the sub bin boundaries
                 std::vector<float>::const_iterator mbvalue = m_boundaries.begin();
                 for (; mbvalue != m_boundaries.end(); ++mbvalue){
                     // should define numerically stable
                     if (fabs((*mbvalue)-sBinMin)<10e-10){
                         // copy the sub bin boundaries into the vector
                         m_totalBoundaries.insert(m_totalBoundaries.begin(),subBinBoundaries.begin(),subBinBoundaries.end());
                         ++mbvalue;
                     } else
                         m_totalBoundaries.push_back(*mbvalue);
                 }
             } else { // (B) multiplicative sub structure
                 // every bin is just repaced by the sub binning structure
                 m_totalBins = m_bins*subBinningData->bins();
                 m_totalBoundaries.reserve(m_totalBins+1);
                 // get the sub bin boundaries if there are any
                 const std::vector<float>& subBinBoundaries = subBinningData->boundaries();
                 // create the boundary vector
                 m_totalBoundaries.push_back(min);
                 for (size_t ib = 0; ib < m_bins; ++ib){
                     float offset = ib*step;
                     for (size_t isb = 1; isb < subBinBoundaries.size(); ++isb)
                         m_totalBoundaries.push_back(offset+subBinBoundaries[isb]);
                 }
             }
             // sort the total boundary vector
             std::sort(m_totalBoundaries.begin(),m_totalBoundaries.end());
         }
      }

     /** Equidistant search : equidist 0 */
     static size_t searchEquidstantWithBoundary(float value, const BinningData& bData) {

       int bin = ((value-bData.min)/bData.step);
       // special treatment of the 0 bin for closed
       if ( bData.option == closed){
           if (value < bData.min) return (bData.m_bins-1);
           if (value > bData.max) return 0;
       }
       // if outside boundary : return boundary for open, opposite bin for closed
       bin = bin < 0 ? ( ( bData.option == open) ? 0 : (bData.m_bins-1) ) : bin ;
       return size_t((bin <= int(bData.m_bins-1)) ? bin : ( (bData.option == open) ? (bData.m_bins-1) : 0 ));
     }

     /** Linear search in vector - superior in O(10) searches: arbitraty 2*/
     static size_t searchInVectorWithBoundary(float value, const BinningData& bData)
     {
       if (bData.binvalue == binPhi) while (value<bData.m_boundaries[0]) value += 2*M_PI;
       if (bData.binvalue == binPhi) while (value>bData.max) value -= 2*M_PI;
       // lower boundary
       if ( value <= bData.m_boundaries[0] ) {
	       return ( bData.option==closed ) ? (bData.m_bins-1) : 0;
       }
       // higher boundary
       if ( value >= bData.max ) return ( bData.option==closed ) ? 0 : (bData.m_bins-1);
       // search
       auto vIter = bData.m_boundaries.begin();
       size_t bin = 0;
       for ( ; vIter !=  bData.m_boundaries.end() ;  ++vIter, ++bin )
	       if ((*vIter) > value) break;
       return (bin-1);
     }

     /** A binary search with underflow/overflow - faster than vector search for O(50) objects*/
     static size_t binarySearchWithBoundary(float value, const BinningData& bData)
     {
       // Binary search in an array of n values to locate value
       //!< @TODO exchange acos(-1) with a const value
       if (bData.binvalue == binPhi) while (value<bData.m_boundaries[0]) value+= 2*acos(-1.);
       if (bData.binvalue == binPhi) while (value>bData.max) value-= 2*acos(-1.);
       // underflow
       if ( value <= bData.m_boundaries[0] ) return ( bData.option==closed ) ? (bData.m_bins-1) : 0;
       size_t nabove, nbelow, middle;
       // overflow
       nabove = bData.m_boundaries.size()+1;
       if ( value >=  bData.max) return ( bData.option==closed ) ? 0 : nabove-2;
       // binary search
       nbelow = 0;
       while(nabove-nbelow > 1) {
	    middle = (nabove+nbelow)/2;
	    if (value == bData.m_boundaries[middle-1]) return middle-1;
	    if (value  < bData.m_boundaries[middle-1]) nabove = middle;
	    else nbelow = middle;
       }
       return nbelow-1;
     }

     /** Search in mixed vector - linear in O-10 m_bins, otherwise binary */
     static size_t searchInVectorWithMixedBoundary(std::pair<float,float> val, const BinningData& bData)
     {
       if ( (bData.hbounds[0].first==0 ? val.first:val.second) < bData.hbounds[0].second ) return ( bData.option==closed ) ? (bData.m_bins-1) : 0;
       if ( (bData.hbounds.back().first==0 ? val.first:val.second)  >= bData.max ) return ( bData.option==closed ) ? 0 : (bData.m_bins-1);

       if ( bData.hbounds.size()<10 ) {
	     auto vBeg = bData.hbounds.begin();
	     auto vIter = vBeg+1;
	     for (  ; vIter !=  bData.hbounds.end() ;  ++vIter )
	         if ((*vIter).second > ((*vIter).first==0 ? val.first:val.second) ) break;
	     return ( vIter!= bData.hbounds.end() ? vIter-vBeg-1 : bData.m_bins-1 ) ;
       }

       // Binary search in an array of n values to locate value
       size_t nabove, nbelow, middle;
       nabove = bData.hbounds.size();
       // binary search
       nbelow = 0;
       while(nabove-nbelow > 1) {
          middle = (nabove+nbelow)/2;
          float valm = bData.hbounds[middle].first==0 ? val.first:val.second;
          if (valm == bData.hbounds[middle].second) { nbelow = middle; break;}
          if (valm  < bData.hbounds[middle].second) nabove = middle;
          else nbelow = middle;
       }

       if (nbelow > bData.m_bins-1 ) return bData.m_bins-1;
       return nbelow;
     }
   };

}

#endif // ACTS_GEOMETRYUTILS_BINNINGDATA_H

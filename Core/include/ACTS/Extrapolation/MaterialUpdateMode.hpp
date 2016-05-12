#ifndef ACTS_EXTRAPOLATIONUTILS_MATERIALUPDATEMODE_H
#define ACTS_EXTRAPOLATIONUTILS_MATERIALUPDATEMODE_H 1

namespace Acts {

   /** This is a steering enum to force the material update 
       it can be:
        (1)  addNoise
        (-1) removeNoise
       Second is mainly for vertex reconstruction, but potentially dangeraous.
     */

   enum MaterialUpdateMode 
   {
          addNoise    =  1,
          removeNoise = -1 
   };
   
}

#endif

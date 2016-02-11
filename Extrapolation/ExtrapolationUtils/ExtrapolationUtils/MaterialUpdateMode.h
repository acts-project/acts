#ifndef ATS_EXTRAPOLATIONUTILS_MATERIALUPDATEMODE_H
#define ATS_EXTRAPOLATIONUTILS_MATERIALUPDATEMODE_H 1

namespace Ats {

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

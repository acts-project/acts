///////////////////////////////////////////////////////////////////
// BinningType.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYUTILS_BINNINGTYPE_H
#define ACTS_GEOMETRYUTILS_BINNINGTYPE_H 1

// STL include(s)
#include <vector>

namespace Acts {

  /** @enum BinningType, BinningOption & BinningAccess

     - BinningType:

       Enumeration to qualify the binning type for the use of the
       LayerArrayCreator and the TrackingVolumeArrayCreator

      - BinningOption:
        open:   [0,max]
        closed:  0 -> nextbin -> max -> 0

      - BinningValue
        necessary access to global positions

      @author Andreas.Salzburger@cern.ch
     */
  enum BinningType { equidistant, arbitrary };

  /** enum BinValue */
  enum BinningOption { open, closed };

  /**  how to take the global / local position */
  enum BinningValue { binX, binY, binZ, binR, binPhi, binRPhi, binH, binEta };

  /** screen output option */
  static std::vector< std::string > binningValueNames = { "binX", "binY", "binZ", "binR", "binPhi", "binRPhi", "binH", "binEta" };

}
#endif // ACTS_GEOMETRYUTILS_BINNINGTYPE_H


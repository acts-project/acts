///////////////////////////////////////////////////////////////////
// TransportJacobian.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef TRKEXUTILS_TRANSPORTJACOBIAN_H
#define TRKEXUTILS_TRANSPORTJACOBIAN_H


// CLHEP
#include "Algebra/AlgebraDefinitions.h"

class MsgStream;

namespace Ats {

/** @class TransportJacobian
  
   This class represents the jacobian for transforming initial track parameters
   to new track parameters after propagation to new surface.
   Initial track parameters:           loc1  ,loc2  ,Phi  ,Theta  ,qp
   Track parameters after propagation: lol1_n,loc2_n,Phi_n,Theta_n,qp_n

   Jacobian is matrix (5x5) with derivatives

          0                1                2               3                 4 
   0 d loc1_n /d loc1 d loc1_n /d loc2 d loc1_n /d Phi d loc1_n /d Theta d loc1_n /d qp
   1 d loc2_n /d loc1 d loc2_n /d loc2 d loc2_n /d Phi d loc2_n /d Theta d loc2_n /d qp
   2 d Phi_n  /d loc1 d Phi_n  /d loc2 d Phi_n  /d Phi d Phi_n  /d Theta d Phi_n  /d qp
   3 d Theta_n/d loc1 d Theta_n/d loc2 d Theta_n/d Phi d Theta_n/d Theta d Theta_n/d qp
   4 d qp_n   /d loc1 d qp_n   /d loc2 d qp_n   /d Phi d qp   _n/d Theta d qp_n   /d qp
  
   ^ ---> second index     
   |
   | first index 


   @author Igor.Gavrilenko@cern.ch 
  
  */

  class TransportJacobian : public AtsMatrixD<5,5>
    {
    public:
    
      /** Constructor */
      TransportJacobian(const double*   );
      TransportJacobian(const AtsMatrixD<5,5>&);
      
      /** Destructor */
      ~TransportJacobian(){};
      
    private:
      
    };
    

  /**Overload of << operator for both, MsgStream and std::ostream for debug output*/
 
  MsgStream& operator << ( MsgStream& sl, const TransportJacobian& jac);
  std::ostream& operator << ( std::ostream& sl, const TransportJacobian& jac); 
  
} // end of namespace

#endif // TRKEXUTILS_TRANSPORTJACOBIAN_H

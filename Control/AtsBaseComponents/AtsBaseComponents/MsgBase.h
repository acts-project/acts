///////////////////////////////////////////////////////////////////
// MsgMacros.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef BASE_COMPINTERFACES_H
#define BASE_COMPINTERFACES_H 1

// STL includes
#include <iosfwd>
#include <string>

// framework includes
#include "GaudiKernel/IMessageSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/Property.h"

namespace Bc {
    
    /** @class MsgBase 
    
         Mixin class to provide easy @c MsgStream access and capabilities.
         One usually inherits from this class and use it like so:
         @code
          void some_method (MsgBase& o)
          { o.msg() << "foo" << endreq; }
         @endcode
    
         @author Sebastian Binet, adapted by Andreas.Salzburger@cern.ch
    
     */
    class MsgBase { 
    
      public: 
        /** Constructor - forbidden */  
        MsgBase() = delete;
        
        /** Copy constructor - forbidden */
        MsgBase( const MsgBase& rhs ) = delete; 
        
        /** Assignmnet operator - forbidden */
        MsgBase& operator=( const MsgBase& rhs ) = delete;
    
        /** Constructor from MsgStream Service */
        MsgBase (IMessageSvc* msgSvc, const std::string& name);
          
        /** Destructor */
        virtual ~MsgBase(); 
          
          
        /** @brief Test the output level
         *  @param lvl The message level to test against
         *  @return boolean Indicting if messages at given level will be printed
         *  @retval true Messages at level "lvl" will be printed
         */
        bool msgLvl (const MSG::Level lvl) const;
          
        /** The standard message stream.
         *  Returns a reference to the default message stream
         *  May not be invoked before sysInitialize() has been invoked.
         */
        MsgStream& msg() const;
          
        /** The standard message stream.
         *  Returns a reference to the default message stream
         *  May not be invoked before sysInitialize() has been invoked.
         */
        MsgStream& msg (const MSG::Level lvl) const;
    
    private: 
    
      /// Default constructor:
     
      ///////////////////////////////////////////////////////////////////
      // Private data:
      ///////////////////////////////////////////////////////////////////
     private: 
    
      /// MsgStream instance (a std::cout like with print-out levels)
      mutable MsgStream m_msg;
    }; 
    
    ///////////////////////////////////////////////////////////////////
    // Inline methods:
    ///////////////////////////////////////////////////////////////////
    //std::ostream& operator<<( std::ostream& out, const MsgBase& o );
    
    inline bool MsgBase::msgLvl (const MSG::Level lvl) const {
     if (m_msg.level() <= lvl) { 
         m_msg << lvl; return true;
     } else {
       return false;
     }
    }
    
    inline MsgStream& MsgBase::msg() const 
    { return m_msg; }
    
    inline MsgStream& MsgBase::msg (const MSG::Level lvl) const 
    { return m_msg << lvl; }

}

#endif 
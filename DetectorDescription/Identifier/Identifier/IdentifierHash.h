/***************************************************************************
 Identifier package
 -----------------------------------------
 Copyright (C) 2002 by ATLAS Collaboration
 ***************************************************************************/

//<doc><file>	$Id: IdentifierHash.h,v 1.4 2003-11-26 12:28:21 schaffer Exp $
//<version>	$Name: not supported by cvs2svn $

#ifndef IDENTIFIER_IDENTIFIERHASH_H
# define IDENTIFIER_IDENTIFIERHASH_H

#define __IDENTIFIER_NOACCESSORS__
#ifdef __IDENTIFIER_NOACCESSORS__
#include "GaudiKernel/MsgStream.h"
#endif

//<<<<<< INCLUDES                                                       >>>>>>
//<<<<<< PUBLIC DEFINES                                                 >>>>>>
//<<<<<< PUBLIC CONSTANTS                                               >>>>>>
//<<<<<< PUBLIC TYPES                                                   >>>>>>
//<<<<<< PUBLIC VARIABLES                                               >>>>>>
//<<<<<< PUBLIC FUNCTIONS                                               >>>>>>
//<<<<<< CLASS DECLARATIONS                                             >>>>>>

/**
 **
 **  ---------------------------------------------------
 **  
 **  IdentifierHash :
 **  	
 **  This is a "hash" representation of an Identifier. This encodes a
 **  32 bit index which can be used to look-up "Identifiable"s stored
 **  in a simple vector. It is intended to be a continuous hash,
 **  i.e. it runs from 0 to N-1, where there are N different possible
 **  values for an Identifier(32) within a specific context.
 **
 **  IdentifierHashes are created by default in an invalid state which
 **  can be checked with "is_valid" method. This allows some error
 **  checking.
 **  
 **  ---------------------------------------------------
 */  
class IdentifierHash
{
public:


    ///----------------------------------------------------------------
    /// Define public typedefs
    ///----------------------------------------------------------------
    typedef unsigned int	value_type;

    ///----------------------------------------------------------------
    /// Constructors
    ///----------------------------------------------------------------

    /// Default constructor
    IdentifierHash ();

    /// Copy constructor
    IdentifierHash (const IdentifierHash& other);

    /// Assignment
    IdentifierHash& operator= (const IdentifierHash& other);

    /// Initialization with value
    IdentifierHash (value_type value);

    ///----------------------------------------------------------------
    /// Accessors
    ///----------------------------------------------------------------

    /// Get the value 
    operator unsigned int	(void) const;
    unsigned int value          (void) const;

    ///----------------------------------------------------------------
    /// Error management
    ///----------------------------------------------------------------

    /// Check if id is in a valid state
    bool is_valid () const;

    ///----------------------------------------------------------------
    /// Modifications
    ///----------------------------------------------------------------

    /// Assignment operators
    IdentifierHash& operator = (value_type value);
    IdentifierHash& operator += (unsigned int value);
    IdentifierHash& operator -= (unsigned int value);

private:

    typedef enum {
        max_value = 0xFFFFFFFF
    } max_value_type;

    //----------------------------------------------------------------
    // The actual identifier data.
    //----------------------------------------------------------------
    value_type m_value;
};
//-----------------------------------------------



//<<<<<< INLINE PUBLIC FUNCTIONS                                        >>>>>>
//<<<<<< INLINE MEMBER FUNCTIONS                                        >>>>>>


//-----------------------------------------------
inline IdentifierHash::IdentifierHash ()
    : m_value(max_value)
{}

//-----------------------------------------------
inline IdentifierHash::IdentifierHash (const IdentifierHash& other)
    : m_value(other.m_value)
{}

//-----------------------------------------------
inline IdentifierHash& IdentifierHash::operator= (const IdentifierHash& other)
{
  if (this != &other)
    m_value = other.m_value;
  return *this;
}

//-----------------------------------------------
inline IdentifierHash::IdentifierHash (value_type value)
    : m_value(value)
{}

//-----------------------------------------------
inline IdentifierHash&
IdentifierHash::operator = (value_type value)
{
    m_value = value;
    return (*this);
}

//-----------------------------------------------
inline IdentifierHash& 				     
IdentifierHash::operator += (unsigned int value)
{
    m_value += value;
    return (*this);
}

//-----------------------------------------------
inline IdentifierHash& 
IdentifierHash::operator -= (unsigned int value)
{
    m_value = (m_value > value) ? m_value - value : 0;
    return (*this);
}

//-----------------------------------------------
inline IdentifierHash::operator unsigned int (void) const
{
    return (m_value);
}

//-----------------------------------------------
inline unsigned int IdentifierHash::value (void) const
{
    return (m_value);
}

//-----------------------------------------------
inline bool 
IdentifierHash::is_valid () const
{
    return (!(max_value == m_value));
}

#ifdef __IDENTIFIER_NOACCESSORS__
inline MsgStream& operator << (MsgStream& f, const IdentifierHash& id)
{
  f << id.value();
  return f;
}

inline std::ostream& operator << (std::ostream& os, const IdentifierHash& id)
{
  os << id.value();
  return os;
}

#endif /* __IDENTIFIER_NOACCESORS__ */

#endif // IDENTIFIER_IDENTIFIERHASH_H

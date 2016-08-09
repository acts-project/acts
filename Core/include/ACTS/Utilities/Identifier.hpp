// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// Identifier.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_CORE_IDENTIFIER_H
#define ACTS_CORE_IDENTIFIER_H 1

#ifdef ACTS_CORE_IDENTIFIER_PLUGIN
#include ACTS_CORE_IDENTIFIER_PLUGIN
#else

#define IDENTIFIER_TYPE unsigned long long
#define IDENTIFIER_DIFF_TYPE long long

#include <string>

/// @class Identifier
///
/// minimum implementation of an Identifier,
/// please use the ACTS_CORE_IDENTIFIER_PLUGING in to use instead if
/// another type of Identifier is needed
///
class Identifier
{
public:
  ///----------------------------------------------------------------
  /// Define public typedefs
  ///----------------------------------------------------------------
  typedef Identifier           id_type;
  typedef IDENTIFIER_TYPE      value_type;
  typedef IDENTIFIER_DIFF_TYPE diff_type;
  typedef IDENTIFIER_TYPE      size_type;

  typedef enum {
    NBITS    = sizeof(value_type) * 8,  // bits per byte
    MAX_BIT  = (static_cast<value_type>(1) << (NBITS - 1)),
    ALL_BITS = ~(static_cast<value_type>(0))
  } bit_defs;

  ///----------------------------------------------------------------
  /// Constructors
  ///----------------------------------------------------------------

  /// Default constructor
  Identifier();

  /// Constructor from value_type
  explicit Identifier(value_type value);

  /// Copy constructor
  Identifier(const Identifier& other);

  ///----------------------------------------------------------------
  /// Modifications
  ///----------------------------------------------------------------
  Identifier&
  operator|=(value_type value);
  Identifier&
  operator&=(value_type value);

  ///----------------------------------------------------------------
  /// Assignment operator
  ///----------------------------------------------------------------
  Identifier&
  operator=(const Identifier& old);
  Identifier&
  operator=(value_type value);

  ///----------------------------------------------------------------
  /// Comparison operators
  ///----------------------------------------------------------------
  bool
  operator==(const Identifier& other) const;
  bool
  operator!=(const Identifier& other) const;
  bool
  operator<(const Identifier& other) const;
  bool
  operator>(const Identifier& other) const;
  bool
  operator<=(const Identifier& other) const;
  bool
  operator>=(const Identifier& other) const;

  /// Check if id is in a valid state
  bool
  is_valid() const;

private:
  //----------------------------------------------------------------
  // The compact identifier data.
  //----------------------------------------------------------------
  value_type m_id;

  typedef enum {
    // max_value = 0xFFFFFFFFFFFFFFFFULL
    max_value = ~(static_cast<value_type>(0))
  } max_value_type;
};
//-----------------------------------------------

//<<<<<< INLINE MEMBER FUNCTIONS                                        >>>>>>

// Constructors
//-----------------------------------------------
inline Identifier::Identifier() : m_id(max_value)
{
}

//-----------------------------------------------
inline Identifier::Identifier(const Identifier& other) : m_id(other.m_id)
{
}

//-----------------------------------------------
inline Identifier::Identifier(value_type value) : m_id(value)
{
}

// Modifications
//-----------------------------------------------

inline Identifier&
Identifier::operator=(const Identifier& other)
{
  if (&other != this) {
    m_id = other.m_id;
  }
  return (*this);
}

inline Identifier&
Identifier::operator=(value_type value)
{
  m_id = value;
  return (*this);
}

inline Identifier&
Identifier::operator|=(value_type value)
{
  m_id |= value;
  return (*this);
}

inline Identifier&
Identifier::operator&=(value_type value)
{
  m_id &= value;
  return (*this);
}

// Comparison operators
//----------------------------------------------------------------
inline bool
Identifier::operator==(const Identifier& other) const
{
  return (m_id == other.m_id);
}

//----------------------------------------------------------------
inline bool
Identifier::operator!=(const Identifier& other) const
{
  return (m_id != other.m_id);
}

//-----------------------------------------------
inline bool
Identifier::operator<(const Identifier& other) const
{
  return (m_id < other.m_id);
}

//-----------------------------------------------
inline bool
Identifier::operator>(const Identifier& other) const
{
  return (m_id > other.m_id);
}

//-----------------------------------------------
inline bool
Identifier::operator<=(const Identifier& other) const
{
  return (m_id <= other.m_id);
}

//-----------------------------------------------
inline bool
Identifier::operator>=(const Identifier& other) const
{
  return (m_id >= other.m_id);
}

//----------------------------------------------------------------
inline bool
Identifier::is_valid() const
{
  return (!(max_value == m_id));
}

#endif  // ACTS_CORE_IDENTIFIER_PLUGIN

#endif  // ACTS_CORE_IDENTIFIER_H

#ifndef __Identifier32_h__
#define __Identifier32_h__

#include <vector>
#include <string>

/**
 **-----------------------------------------------
 **
 **  Identifier32 is a simple type-safe 32 bit unsigned integer. An
 **  Identifier32 relies on other classes - IdHelpers - to encode and
 **  decode its information.
 **  
 **  The default constructor created an Identifier32 an invalid state
 **  which can be check with the "is_valid" method to allow some error
 **  checking.
 **  
 **-----------------------------------------------
 */
class Identifier32
{
public:


    ///----------------------------------------------------------------
    /// Define public typedefs
    ///----------------------------------------------------------------
    typedef Identifier32                id_type;
    typedef unsigned int                value_type;
    typedef unsigned int                size_type;

    ///----------------------------------------------------------------
    /// Constructors
    ///----------------------------------------------------------------

    /// Default constructor
    Identifier32 ();

    /// Constructor from value_type
    explicit Identifier32 (value_type value);

    /// Copy constructor
    Identifier32 (const Identifier32& other);

    /// Assignment.
    Identifier32& operator= (const Identifier32& other);

    ///----------------------------------------------------------------
    /// Modifications
    ///----------------------------------------------------------------

    /// Assignment operator
    Identifier32& operator = (value_type value);

    /// Bitwise operations 
    Identifier32& operator |= (value_type value);
    Identifier32& operator &= (value_type value);

    /// Reset to invalid state
    void clear ();

    ///----------------------------------------------------------------
    /// Accessors
    ///----------------------------------------------------------------

    /// Get the compact id
    value_type  get_compact  (void) const;

    ///----------------------------------------------------------------
    /// Comparison operators
    ///----------------------------------------------------------------
  
    bool operator ==    (const Identifier32& other) const;
    bool operator !=    (const Identifier32& other) const;
    bool operator <     (const Identifier32& other) const;
    bool operator >     (const Identifier32& other) const;
    bool operator <=    (const Identifier32& other) const;
    bool operator >=    (const Identifier32& other) const;

    ///----------------------------------------------------------------
    /// Error management
    ///----------------------------------------------------------------

    /// Check if id is in a valid state
    bool is_valid () const;

    ///----------------------------------------------------------------
    /// Utilities
    ///----------------------------------------------------------------

    /// Provide a string form of the identifier - hexadecimal
    std::string  getString() const;

    /// Print out in hex form
    void show () const;

private:

    typedef enum {
        max_value = 0xFFFFFFFF
    } max_value_type;

    //----------------------------------------------------------------
    // The compact identifier data.
    //----------------------------------------------------------------
    value_type m_id;

};
//-----------------------------------------------










//<<<<<< INLINE MEMBER FUNCTIONS                                        >>>>>>


// Constructors
//-----------------------------------------------
inline Identifier32::Identifier32 ()
    : m_id(max_value)
{}

//-----------------------------------------------
inline Identifier32::Identifier32 (const Identifier32& other)
    : m_id(other.m_id)
{}

//-----------------------------------------------
inline Identifier32& Identifier32::operator= (const Identifier32& other)
{
  if (this != &other)
    m_id = other.m_id;
  return *this;
}

//-----------------------------------------------
inline Identifier32::Identifier32 (value_type value)
    : m_id(value)
{}

// Modifications
//-----------------------------------------------

inline Identifier32&
Identifier32::operator = (value_type value)
{
    m_id = value;
    return (*this);
}

inline Identifier32&                                   
Identifier32::operator |= (unsigned int value)
{
    m_id |= value;
    return (*this);
}

inline Identifier32& 
Identifier32::operator &= (unsigned int value)
{
    m_id &= value;
    return (*this);
}

inline void 
Identifier32::clear () 
{
    m_id = max_value;
}

// Accessors

inline Identifier32::value_type  Identifier32::get_compact  (void) const
{
    return (m_id);
}

// Comparison operators
//----------------------------------------------------------------
inline bool 
Identifier32::operator == (const Identifier32& other) const
{
    return (m_id == other.m_id);
}

//----------------------------------------------------------------
inline bool 
Identifier32::operator != (const Identifier32& other) const
{
    return (m_id != other.m_id);
}

//-----------------------------------------------
inline bool 
Identifier32::operator < (const Identifier32& other) const
{
    return (m_id < other.m_id);
}

//-----------------------------------------------
inline bool 
Identifier32::operator > (const Identifier32& other) const
{
    return (m_id > other.m_id);
}

//-----------------------------------------------
inline bool 
Identifier32::operator <= (const Identifier32& other) const
{
    return (m_id <= other.m_id);
}

//-----------------------------------------------
inline bool 
Identifier32::operator >= (const Identifier32& other) const
{
    return (m_id >= other.m_id);
}

inline bool 
Identifier32::is_valid () const
{
    return (!(max_value == m_id));
}

#endif

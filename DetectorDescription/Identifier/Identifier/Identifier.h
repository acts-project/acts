#ifndef __Identifier_h__
#define __Identifier_h__

// #define __IDENTIFIER_64BIT__
#define __IDENTIFIER_NOACCESSORS__

#ifndef __IDENTIFIER_64BIT__
#define IDENTIFIER_TYPE unsigned int
#define IDENTIFIER_DIFF_TYPE int
#define IDENTIFIER_PCODE ""
#else /* __IDENTIFIER_64BIT__ */
#define IDENTIFIER_TYPE unsigned long long
#define IDENTIFIER_DIFF_TYPE long long
#define IDENTIFIER_PCODE "ll"
#endif /* __IDENTIFIER_64BIT__ */

#ifdef __IDENTIFIER_NOACCESSORS__
#include "GaudiKernel/MsgStream.h"
#endif

#include "Identifier/Identifier32.h"

#include <iostream>
#include <boost/io/ios_state.hpp>
#include <vector>
#include <string>

/**
 **-----------------------------------------------
 **
 **  Identifier is a simple type-safe 64 bit unsigned integer. An
 **  Identifier relies on other classes - IdHelpers - to encode and
 **  decode its information.
 **  
 **  The default constructor created an Identifier an invalid state
 **  which can be check with the "is_valid" method to allow some error
 **  checking.
 **  
 **-----------------------------------------------
 */
class Identifier
{
public:


    ///----------------------------------------------------------------
    /// Define public typedefs
    ///----------------------------------------------------------------
    typedef Identifier                  id_type;
    typedef IDENTIFIER_TYPE             value_type;
    typedef IDENTIFIER_DIFF_TYPE        diff_type;
    typedef IDENTIFIER_TYPE             size_type;

    typedef enum
    {
        NBITS = sizeof(value_type) * 8, // bits per byte
        MAX_BIT = (static_cast<value_type>(1) << (NBITS - 1)),
        ALL_BITS = ~(static_cast<value_type>(0))
    } bit_defs;

    ///----------------------------------------------------------------
    /// Constructors
    ///----------------------------------------------------------------

    /// Default constructor
    Identifier ();

    /// Constructor from value_type
    explicit Identifier (value_type value);

    /// Constructor from Identifier32
    Identifier (const Identifier32& other);

#ifdef __IDENTIFIER_64BIT__
    /// Constructor from 32-bit value_type and int
    /// (to avoid common implicit conversions)
    explicit Identifier (Identifier32::value_type value);
    explicit Identifier (int value);
#endif

    /// Copy constructor
    Identifier (const Identifier& other);

    ///----------------------------------------------------------------
    /// Modifications
    ///----------------------------------------------------------------

    /// Assignment operator
    Identifier& operator = (const Identifier& old);
    Identifier& operator = (const Identifier32& old);
    Identifier& operator = (value_type value);
#ifdef __IDENTIFIER_64BIT__
    /// Assignment to avoid common implicit conversions and shift properly
    Identifier& operator = (Identifier32::value_type value);
    Identifier& operator = (int value);
#endif

    /// Bitwise operations 
//#ifndef __IDENTIFIER_64BIT__
private:
    Identifier& operator |= (value_type value);
    Identifier& operator &= (value_type value);
public:
//#endif
    //Identifier& operator |= (const Identifier& other);
    //Identifier& operator &= (const Identifier& other);
    // these methods make sure the proper shift is done
    // (perhaps don't need it - Identifier(Identifier32) is implicit
    //Identifier& operator |= (const Identifier32& other);
    //Identifier& operator &= (const Identifier32& other);

    /// build from a string form - hexadecimal
    void set (const std::string& id);

    /// Reset to invalid state
    void clear ();

    /// Set literal value
    Identifier& set_literal(value_type value);

    ///----------------------------------------------------------------
    /// Accessors
    ///----------------------------------------------------------------

#ifndef __IDENTIFIER_NOACCESSORS__ // get rid of conversion accessors
    /// Get the value 
    operator     value_type      (void) const;
#endif /* __IDENTIFIER_NOACCESSORS__ */

    /// Get the 32-bit version Identifier, will be invalid if >32 bits
    /// needed
    Identifier32 get_identifier32  (void) const;

    /// Get the compact id
    value_type   get_compact  (void) const;

    /// A get_compact functional for use in STL algorithms
    class get_compact_func
    {
    public:
        value_type operator() (const Identifier& id)
            {
                return id.get_compact();
            }
    };

    ///----------------------------------------------------------------
    /// Comparison operators
    ///----------------------------------------------------------------
  
    bool operator ==    (const Identifier& other) const;
    bool operator !=    (const Identifier& other) const;
    bool operator <     (const Identifier& other) const;
    bool operator >     (const Identifier& other) const;
    bool operator <=    (const Identifier& other) const;
    bool operator >=    (const Identifier& other) const;

    //bool operator ==    (const Identifier32& other) const;
    //bool operator !=    (const Identifier32& other) const;
    //bool operator <     (const Identifier32& other) const;
    //bool operator >     (const Identifier32& other) const;
    //bool operator <=    (const Identifier32& other) const;
    //bool operator >=    (const Identifier32& other) const;

#ifdef __IDENTIFIER_NOACCESSORS__
    /// Comparison operators with value_type.
    /// This is a hack, only because GeoAdaptors/GeoMuonHits wants to
    /// to compare explicitly with 0 as a test of whether the identifier
    /// has been constructed properly.  But is_valid() here compares
    /// with max_value, not 0, since presumably it is possible to have
    /// a zero value - just not in muons.
    bool operator ==    (value_type other) const;
    bool operator !=    (value_type other) const;
#ifdef __IDENTIFIER_64BIT__
    bool operator ==    (Identifier32::value_type other) const;
    bool operator ==    (int other) const;
    bool operator !=    (Identifier32::value_type other) const;
    bool operator !=    (int other) const;
#endif /* __IDENTIFIER_64BIT__ */
#endif /* __IDENTIFIER_NOACCESORS__ */

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

    /// extract field from identifier (shift first, then mask)
    value_type extract(size_type shift, size_type mask) const;

    /// extract field(s) by masking first, then shifting
    value_type mask_shift(value_type mask, size_type shift) const;

    /// extract field, no mask
    value_type extract(size_type shift) const;

    // allow IdDict access to the following private methods
    friend class IdDictDictionary;
    friend class IdDictFieldImplementation;
    friend class AtlasDetectorID;
    friend class PixelID;
    
    typedef enum {
        //max_value = 0xFFFFFFFFFFFFFFFFULL
        max_value = ~(static_cast<value_type>(0))
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
inline Identifier::Identifier ()
        : m_id(max_value)
{}

//-----------------------------------------------
inline Identifier::Identifier (const Identifier& other)
        : m_id(other.m_id)
{}

/// Constructor from Identifier32
//-----------------------------------------------
inline Identifier::Identifier (const Identifier32& other)
        : m_id(max_value)
{
    //std::cout << "Identifier(Identifier32) " << other.get_compact() << std::endl;
    if (other.is_valid()) {
#ifndef __IDENTIFIER_64BIT__
        m_id = other.get_compact();
#else /* __IDENTIFIER_64BIT__ */
        m_id = (static_cast<value_type>(other.get_compact()) << 32);
#endif /* __IDENTIFIER_64BIT__ */
    }
}

#ifdef __IDENTIFIER_64BIT__
/// Constructor from Identifier32 value_type (unsigned int)
/// (only use in id64 case since otherwise redundant)
//-----------------------------------------------
inline Identifier::Identifier (Identifier32::value_type value)
        : m_id(max_value)
{
    //std::cout << "Identifier(Identifier32::value_type) " << value << std::endl;
    m_id = (static_cast<value_type>(value) << 32);
}
inline Identifier::Identifier (int value)
        : m_id(max_value)
{
    //std::cout << "Identifier(int) " << value << std::endl;
    m_id = (static_cast<value_type>(value) << 32);
}
#endif /* __IDENTIFIER_64BIT__ */

//-----------------------------------------------
inline Identifier::Identifier (value_type value)
        : m_id(value)
{
    //std::cout << "Identifier(value_type) " << value << std::endl;
#ifdef __IDENTIFIER_64BIT__
    // Print out warning for potential call with value for a 32-bit id
    // I.e. if lower bits are set and no upper bit set
    const value_type upper = 0XFFFFFFFF00000000LL;
    const value_type lower = 0X00000000FFFFFFFFLL;
    const value_type testUpper = value & upper;
    const value_type testLower = value & lower;
    if ( testUpper == 0 && testLower > 0) {
        boost::io::ios_flags_saver ifs(std::cout);
        std::cout << "Identifier::Identifier - WARNING Constructing 64-bit id with empty upper and non-empty lower: " << std::hex << testUpper << " " << testLower << std::endl;
        m_id = (value << 32);
    }
#endif /* __IDENTIFIER_64BIT__ */

}

// Modifications
//-----------------------------------------------

inline Identifier&
Identifier::operator = (const Identifier& old) {
  if (&old != this) {
    //std::cout << "operator=(Identifier) " << old.get_compact() << std::endl;
    m_id = old.m_id;
  }
  return (*this);
}

inline Identifier&
Identifier::operator = (const Identifier32& old) {
    //std::cout << "operator=(Identifier32) " << old.get_compact() << std::endl;
#ifdef __IDENTIFIER_64BIT__
    m_id = (static_cast<value_type>(old.get_compact()) << 32);
#else
    m_id = old.get_compact();
#endif
    return (*this);
}

inline Identifier&
Identifier::operator = (value_type value)
{
    //std::cout << "operator=(value_type) " << value << std::endl;
#ifdef __IDENTIFIER_64BIT__
    // Print out warning for potential call with value for a 32-bit id
    // I.e. if lower bits are set and no upper bit set
    const value_type upper = 0XFFFFFFFF00000000LL;
    const value_type lower = 0X00000000FFFFFFFFLL;
    const value_type testUpper = value & upper;
    const value_type testLower = value & lower;
    if ( testUpper == 0 && testLower > 0) {
        boost::io::ios_flags_saver ifs(std::cout);
        std::cout << "Identifier::opertor = - WARNING Constructing 64-bit id with empty upper and non-empty lower: " << std::hex << testUpper << " " << testLower << std::endl;
        m_id = (value << 32);
        return (*this);
    }
#endif /* __IDENTIFIER_64BIT__ */

    m_id = value;
    return (*this);
}

#ifdef __IDENTIFIER_64BIT__
inline Identifier&
Identifier::operator = (Identifier32::value_type value)
{
    //std::cout << "operator=(Identifier32::value_type) " << value << std::endl;
    m_id = static_cast<value_type>(value) << 32;
    return (*this);
}
inline Identifier&
Identifier::operator = (int value)
{
    //std::cout << "operator=(int) " << value << std::endl;
    m_id = static_cast<value_type>(value) << 32;
    return (*this);
}
#endif /* __IDENTIFIER_64BIT__ */

//#ifndef __IDENTIFIER_64BIT__
inline Identifier&                                   
Identifier::operator |= (value_type value)
{
    m_id |= value;
    return (*this);
}

inline Identifier& 
Identifier::operator &= (value_type value)
{
    m_id &= value;
    return (*this);
}
//#endif /* __IDENTIFIER_64BIT__ */

//inline Identifier&                                   
//Identifier::operator |= (const Identifier& other)
//{
//    m_id |= other.get_compact();
//    return (*this);
//}
//
//inline Identifier& 
//Identifier::operator &= (const Identifier& other)
//{
//    m_id &= other.get_compact();
//    return (*this);
//}

//inline Identifier&                                   
//Identifier::operator |= (const Identifier32& other)
//{
//    m_id |= (static_cast<value_type>(other.get_compact()) << 32);
//    return (*this);
//}
//
//inline Identifier& 
//Identifier::operator &= (const Identifier32& other)
//{
//    m_id &= (static_cast<value_type>(other.get_compact()) << 32);
//    return (*this);
//}

inline Identifier&
Identifier::set_literal (value_type value)
{
    m_id = value;
    return (*this);
}

inline void 
Identifier::clear () 
{
    m_id = max_value;
}

inline Identifier::value_type Identifier::extract(
    Identifier::size_type shift, Identifier::size_type mask) const {
    return (m_id >> shift) & static_cast<Identifier::value_type>(mask);
}

inline Identifier::value_type Identifier::mask_shift(
    Identifier::value_type mask, Identifier::size_type shift) const {
    return (m_id & mask) >> shift;
}

inline Identifier::value_type Identifier::extract(
    Identifier::size_type shift) const {
    return (m_id >> shift);
}

#ifndef __IDENTIFIER_NOACCESSORS__
// Accessors

inline Identifier::operator Identifier::value_type (void) const
{
    return (m_id);
}
#endif /* __IDENTIFIER_NOACCESSORS__ */
                                             
inline Identifier32 Identifier::get_identifier32  (void) const
{
#ifndef __IDENTIFIER_64BIT__
    return (Identifier32(m_id));
#else /* __IDENTIFIER_64BIT__ */
    // test for bit set in lower 32
    if (extract(0,0xFFFFFFFF)) return (Identifier32());
    return (Identifier32(extract(32)));
#endif /* __IDENTIFIER_64BIT__ */
}


inline Identifier::value_type  Identifier::get_compact  (void) const
{
    return (m_id);
}

// Comparison operators
//----------------------------------------------------------------
inline bool 
Identifier::operator == (const Identifier& other) const
{
    return (m_id == other.m_id);
}

//----------------------------------------------------------------
inline bool 
Identifier::operator != (const Identifier& other) const
{
    return (m_id != other.m_id);
}

//-----------------------------------------------
inline bool 
Identifier::operator < (const Identifier& other) const
{
    return (m_id < other.m_id);
}

//-----------------------------------------------
inline bool 
Identifier::operator > (const Identifier& other) const
{
    return (m_id > other.m_id);
}

//-----------------------------------------------
inline bool 
Identifier::operator <= (const Identifier& other) const
{
    return (m_id <= other.m_id);
}

//-----------------------------------------------
inline bool 
Identifier::operator >= (const Identifier& other) const
{
    return (m_id >= other.m_id);
}

//----------------------------------------------------------------
//inline bool 
//Identifier::operator == (const Identifier32& other) const
//{
//    return (this == Identifier(other));
//}
//
//----------------------------------------------------------------
//inline bool 
//Identifier::operator != (const Identifier32& other) const
//{
//    return (this != Identifier(other));
//}
//
//-----------------------------------------------
//inline bool 
//Identifier::operator < (const Identifier32& other) const
//{
//    return (this < Identifier(other));
//}
//
//-----------------------------------------------
//inline bool 
//Identifier::operator > (const Identifier32& other) const
//{
//    return (this > Identifier(other));
//}
//
//-----------------------------------------------
//inline bool 
//Identifier::operator <= (const Identifier32& other) const
//{
//    return (this <= Identifier(other));
//}
//
//-----------------------------------------------
//inline bool 
//Identifier::operator >= (const Identifier32& other) const
//{
//    return (this >= Identifier(other));
//}

#ifdef __IDENTIFIER_NOACCESSORS__
//----------------------------------------------------------------
inline bool 
Identifier::operator == (Identifier::value_type other) const
{
    return (m_id == other);
}

inline bool 
Identifier::operator != (Identifier::value_type other) const
{
    return (m_id != other);
}

#ifdef __IDENTIFIER_64BIT__
inline bool 
Identifier::operator == (Identifier32::value_type other) const
{
    return ((*this) == Identifier(other));
}

inline bool 
Identifier::operator == (int other) const
{
    return ((*this) == Identifier(other));
}

inline bool 
Identifier::operator != (Identifier32::value_type other) const
{
    return ((*this) != Identifier(other));
}

inline bool 
Identifier::operator != (int other) const
{
    return ((*this) != Identifier(other));
}

#endif /* __IDENTIFIER_64BIT__ */

/// This is for logging

inline MsgStream& operator << (MsgStream& f, const Identifier& id)
{
    f << id.getString();
    return f;
}

inline std::ostream& operator << (std::ostream& os, const Identifier& id)
{
    os << id.getString();
    return os;
}

#endif /* __IDENTIFIER_NOACCESORS__ */

inline bool 
Identifier::is_valid () const
{
    return (!(max_value == m_id));
}

#endif

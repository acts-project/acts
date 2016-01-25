#ifndef __ExpandedIdentifier_h__
#define __ExpandedIdentifier_h__

#include <vector>
#include <string>

//-----------------------------------------------
//
//
//  
//  ExpandedIdentifier :
//  
//  Stores a set of numbers.
//  
//  ---------------------------------------------------
//  
//  ------------------------
//  Possible operations :
//  ------------------------
//
//  ------------------------
//  Constructors:
//  ------------------------
//
//  ExpandedIdentifier ()                         : default constructor
//
//  ExpandedIdentifier (const std::string& text)  : build an identifier from a textual
//                                                  description following the syntax 
//        <value> [ "/" <value> ... ]
//
//  ExpandedIdentifier (const ExpandedIdentifier& other)  : copy constructor
//
//  ExpandedIdentifier (const ExpandedIdentifier& other, 
//                      size_type start)          : subset constructor
//
//  ------------------------
//  Initialisations:
//  ------------------------
//
//  void clear ()                         : clears up the identifier
//
//  void set (const std::string& text)    : set from a textual description
//
//  ------------------------
//  Modifications:
//  ------------------------
//
//  void add (element_type value)         : appends a numeric value to
//                                          an identifier (adds a field).
//
//  ExpandedIdentifier& operator << (element_type value)   : appends a numeric value to
//                                          an identifier (adds a field).
//
//  ------------------------
//  Accessors:
//  ------------------------
//
//  size_type fields ()                   : returns the number of fields 
//                                          currently stored into the 
//                                          identifier.
//
//  element_type operator [] (size_type field) : gets the value stored at the
//                                          specified field number.
//
//  ------------------------
//  Comparison operators:
//  ------------------------
//
//  operator < (id_type& other)          : absolute comparison, 
//
//                                         e.g. :
//
//                                       /1/2 < /1/3
//                                       /1   < /1/3
//                                       /1 is implicitly /1/0...
//  
//  prefix_less (id_type& other)         : comparison on the equal length 
//                                         prefix parts of two ids, 
//
//                                         e.g. :
//
//                                       /1/2 < /1/3,
//                                          but now 
//                                       /1 == /1/3  
//
//
//  error_code last_error ()             : returns the last error code 
//                                         produced by the most recent 
//                                         identifier operation. The possible
//                                         values are :
//
//                                             none
//                                             bad_parameter
//                                             field_not_found
//
//  const std::string last_error_text () : returns a text describing the 
//                                         last error.
//
//  ----------------------------------------------------
//
//  Example of how to use an identifier :
//
//  #include <ExpandedIdentifier.h>
//
//  ExpandedIdentifier id;
//
//  id << 125 << 236 << 306 << 2222;
//  if (id.last_error () != ExpandedIdentifier::none)
//    {
//      cout << "Error : " << id.last_error_text () << endl;
//    }
//
//  for (size_type i = 0; i < id.fields (); ++i)
//    {
//      cout << "id[" << i << "] = " << id[i] << endl;
//    }
//
//-----------------------------------------------
class ExpandedIdentifier
{
public:


    //----------------------------------------------------------------
    // Define public typedefs
    //----------------------------------------------------------------
  typedef ExpandedIdentifier 		id_type;
  typedef int  				element_type;
  typedef std::vector<element_type> 	element_vector;
  typedef std::vector<element_type>::size_type 	size_type;

  typedef enum
    {
      max_value = 0x3FFFFFFF
    } max_value_type;

    //----------------------------------------------------------------
    // Possible errors that may occur in operations.
    //----------------------------------------------------------------
  typedef enum
    {
      none,                 // success
      bad_parameter,        // bad value for any method parameter
      field_not_found,      // the field number is not found
      
      errors,
      get
    } error_code;

    //----------------------------------------------------------------
    // Constructors
    //----------------------------------------------------------------

    //----------------------------------------------------------------
    // Default constructor
    //----------------------------------------------------------------
  ExpandedIdentifier ();

    //----------------------------------------------------------------
    // Copy constructor and assignment.
    //----------------------------------------------------------------
  ExpandedIdentifier (const ExpandedIdentifier& other);
  ExpandedIdentifier& operator= (const ExpandedIdentifier& other);

    //----------------------------------------------------------------
    // Constructor from a subset of another ExpandedIdentifier
    //----------------------------------------------------------------
  ExpandedIdentifier (const ExpandedIdentifier& other, size_type start);

    //----------------------------------------------------------------
    // Constructor from a textual description
    //----------------------------------------------------------------
  ExpandedIdentifier (const std::string& text);

    //----------------------------------------------------------------
    // Modifications
    //----------------------------------------------------------------

    //----------------------------------------------------------------
    // Append a value into a new field.
    //----------------------------------------------------------------
  void add (element_type value);
  ExpandedIdentifier& operator << (element_type value);
  element_type& operator [] (size_type index);

    //----------------------------------------------------------------
    // build from a textual description
    //----------------------------------------------------------------
  void set (const std::string& text);

    //----------------------------------------------------------------
    // Erase all fields.
    // All previously stored data is lost.
    //----------------------------------------------------------------
  void clear ();

    //----------------------------------------------------------------
    // Accessors
    //----------------------------------------------------------------

    //----------------------------------------------------------------
    // Get the value stored into the specified field.
    //----------------------------------------------------------------
  element_type operator [] (size_type index) const;

    //----------------------------------------------------------------
    // Count the number of fields.
    //----------------------------------------------------------------
  size_type fields () const;

    //----------------------------------------------------------------
    // Comparison operators
    //----------------------------------------------------------------
  
  int operator == (const ExpandedIdentifier& other) const;
  int operator != (const ExpandedIdentifier& other) const;
  int operator < (const ExpandedIdentifier& other) const;
  int operator > (const ExpandedIdentifier& other) const;
  int prefix_less (const ExpandedIdentifier& other) const;

  /**
   *    Test if the shorter of two ids is identical
   *    to the equivalent sub-id extracted from the longer
   */
  int match (const ExpandedIdentifier& other) const;

    //----------------------------------------------------------------
    // Error management
    //----------------------------------------------------------------

    //----------------------------------------------------------------
    // Return the error produced in the last operation
    // The value identifier::none represents the successful condition.
    //----------------------------------------------------------------
  error_code last_error () const;

    //----------------------------------------------------------------
    // Return a textual equivalent to the error produced
    // in the last operation
    //----------------------------------------------------------------
  const std::string last_error_text () const;

    //----------------------------------------------------------------
    // Utilities
    //----------------------------------------------------------------

    //----------------------------------------------------------------
    // Send a textual representation of the identifier using the input format
    //----------------------------------------------------------------
  operator std::string () const;

  void show () const;

private:

    //----------------------------------------------------------------
    // The actual identifier data.
    //----------------------------------------------------------------
  element_vector m_fields;

    //----------------------------------------------------------------
    // Maintains the last error code (shared by all objects).
    //----------------------------------------------------------------
  error_code set_last_error (error_code code = get) const;
};
//-----------------------------------------------

// inline definitions


  // Constructors
//-----------------------------------------------
inline
ExpandedIdentifier::ExpandedIdentifier ()
{
  set_last_error (none);
}

//-----------------------------------------------
inline
ExpandedIdentifier::ExpandedIdentifier (const ExpandedIdentifier& other)
{
  set_last_error (none);

  m_fields = other.m_fields;
}

//-----------------------------------------------
inline
ExpandedIdentifier&
ExpandedIdentifier::operator= (const ExpandedIdentifier& other)
{
  if (this != &other) {
    set_last_error (none);
    m_fields = other.m_fields;
  }
  return *this;
}

//-----------------------------------------------
inline
ExpandedIdentifier::ExpandedIdentifier (const ExpandedIdentifier& other, size_type start)
{
  set_last_error (none);

  if (start < other.fields ())
    {
      element_vector::const_iterator it = other.m_fields.begin ();
      it += start;

      m_fields.insert (m_fields.end (),
                       it, 
                       other.m_fields.end ());
    }
}

  // Modifications
//-----------------------------------------------
inline
void ExpandedIdentifier::add (element_type value)
{
  set_last_error (none);

  // Max size of id levels should be < 10
  if(m_fields.capacity() < 10) m_fields.reserve(10);
  m_fields.push_back (value);
}

//-----------------------------------------------
inline
ExpandedIdentifier& ExpandedIdentifier::operator << (element_type value)
{
  set_last_error (none);

  // Max size of id levels should be < 10
  if(m_fields.capacity() < 10) m_fields.reserve(10);
  m_fields.push_back (value);

  return (*this);
}

//-----------------------------------------------
inline
ExpandedIdentifier::element_type & ExpandedIdentifier::operator [] (size_type index)
{
  set_last_error (none);

  if (index >= m_fields.size ())
    {
      static element_type dummy = 0;

      set_last_error (field_not_found);

      return (dummy);
    }

  return (m_fields[index]);
}

//-----------------------------------------------
inline
void ExpandedIdentifier::clear ()
{
  m_fields.clear ();
}



  // Accessors
//-----------------------------------------------
inline
ExpandedIdentifier::element_type ExpandedIdentifier::operator [] (size_type index) const
{
  set_last_error (none);

  if (index >= m_fields.size ())
    {
      set_last_error (field_not_found);

      return (0);
    }

  return (m_fields[index]);
}

//-----------------------------------------------
inline
ExpandedIdentifier::size_type ExpandedIdentifier::fields () const
{
  return (m_fields.size ());
}

  // Comparison operators

//----------------------------------------------------------------
inline
int ExpandedIdentifier::operator == (const ExpandedIdentifier& other) const
{
  const ExpandedIdentifier& me = *this;
  const size_type my_fields = fields ();
  const size_type other_fields = other.fields ();
  
  if (my_fields != other_fields) return (0);
  
  size_type field = 0;
  for (; field < my_fields; ++field) 
    {
      if (me[field] != other[field]) return (0);
    }

  return (1);
}

//----------------------------------------------------------------
inline
int ExpandedIdentifier::operator != (const ExpandedIdentifier& other) const
{
  const ExpandedIdentifier& me = *this;

  return (!(me == other));
}

//-----------------------------------------------
inline
int ExpandedIdentifier::operator < (const ExpandedIdentifier& other) const
{
  const ExpandedIdentifier& me = *this;
  const size_type my_fields = fields ();
  const size_type other_fields = other.fields ();

  size_type field = 0;
  for (;;)
    {
      if ((field == my_fields) ||
          (field == other_fields))
        {
            // Someone has run out of fields. And up to now my_id ==
            // other_id. If the lengths are different, the following
            // then defines the "shorter" one to be "less than". If
            // the lengths are the same, then the two are NOT "less
            // than".
          return (my_fields < other_fields);
        }

      element_type my_field = me[field];
      element_type other_field = other[field];

      if (my_field < other_field) return (1);
      if (my_field > other_field) return (0);

      field++;
    }

  return (0);
}

//-----------------------------------------------
inline
int ExpandedIdentifier::operator > (const ExpandedIdentifier& other) const
{
  const ExpandedIdentifier& me = *this;

  return (other < me);
}

//----------------------------------------------------------------
inline
int ExpandedIdentifier::match (const ExpandedIdentifier& other) const
{
  const ExpandedIdentifier& me = *this;
  const size_type my_fields = fields ();
  const size_type other_fields = other.fields ();

  const size_type fs = (my_fields < other_fields) ? my_fields : other_fields;
  
  for (size_type field = 0; field < fs; ++field) 
    {
      if (me[field] != other[field]) return (0);
    }

  return (1);
}

  // Error management
//-----------------------------------------------
inline
ExpandedIdentifier::error_code ExpandedIdentifier::last_error () const
{
  return (set_last_error (get));
}

//----------------------------------------------------------------
inline
ExpandedIdentifier::error_code ExpandedIdentifier::set_last_error (error_code code) const
{
  static error_code last = none;
  
  if (code != get)
    {
      last = code;
    }
  
  return (last);
}



#endif

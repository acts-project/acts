
#include "Identifier/ExpandedIdentifier.h"
#include <stdarg.h>
#include <stdio.h>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <iomanip>

//-----------------------------------------------
static void show_vector (const ExpandedIdentifier::element_vector& v)
{
  ExpandedIdentifier::element_vector::const_iterator it;
  bool first = true;

  std::cout << "[";
  for (it = v.begin (); it != v.end (); ++it)
    {
      if (first) first = false;
      else std::cout << ".";

      ExpandedIdentifier::element_type value = *it;
      std::cout << value;

        // if (value >= 0) std::cout << value;
        // else std::cout << "*";
    }
  std::cout << "]";
}


//-----------------------------------------------
ExpandedIdentifier::ExpandedIdentifier (const std::string& text)
{
  set (text);
}


//-----------------------------------------------
void ExpandedIdentifier::set (const std::string& text)
{
  set_last_error (none);

  clear ();
  if (text.size () == 0) return;
  const char* ctext = text.c_str ();

  for (;;)
    {

      const char* sep;

      sep = strchr (ctext, '/');

      int value = 0;
      sscanf (ctext, "%d", &value);

      add ((element_type) value);

      if (sep == 0) break;
      
      ctext = sep + 1;
    }
}


//-----------------------------------------------
int ExpandedIdentifier::prefix_less (const ExpandedIdentifier& other) const
{
  const ExpandedIdentifier& me = *this;
  const size_type my_fields = fields ();
  const size_type other_fields = other.fields ();
  
    //
    // Scan fields up to the less common field number
    //
  size_type field = 0;
  while ((field < my_fields) && 
         (field < other_fields))
    {
      element_type my_field = me[field];
      element_type other_field = other[field];
      
      if (my_field < other_field) return (1);
      if (my_field > other_field) return (0);
      
      field++;
    }
  
  return (0);
}

  // Error management

//-----------------------------------------------
const std::string ExpandedIdentifier::last_error_text () const
{
  static const std::string text[] = {
    "none",
    "bad parameter",
    "field not found"
  };
  
  error_code code = last_error ();
  
  return (text[code]);
}

//-----------------------------------------------
ExpandedIdentifier::operator std::string () const
{
  std::string result;
  char temp[20];

  size_type my_fields = m_fields.size ();

  if (my_fields == 0) return (result);

    // print fields one by one.
  for (size_type field_number = 0; field_number < my_fields; field_number++)
    {
      element_type value = m_fields[field_number];

      if (field_number > 0) result += "/";

      sprintf (temp, "%d", value);
      result += temp;
    }

  return (result);
}

//-----------------------------------------------
void ExpandedIdentifier::show () const
{
  show_vector (m_fields);
}



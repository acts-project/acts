 
 
#include "Identifier/Range.h" 
 
#include <stdio.h> 
#include <string> 
#include <vector> 
#include <algorithm> 
 
#include <limits>
#include <iostream> 
#include <iomanip> 
#include <set>

#include <assert.h> 
 
#ifdef WIN32 
namespace std 
{ 
  template <class T> static const T& min (const T& t1, const T& t2) 
  { 
    return ((t1 < t2) ? t1 : t2); 
  } 
   
  template <class T> static const T& max (const T& t1, const T& t2) 
  { 
    return ((t1 > t2) ? t1 : t2); 
  } 
}; 
#endif 
 
class RangeParser 
{ 
public: 
  typedef Range::size_type size_type; 
 
  RangeParser () 
      { 
      } 
 
  bool run (Range& range, const std::string& text) 
      { 
        range.clear (); 
        size_type pos = 0; 
        return (parse (range, text, pos)); 
      } 
 
private: 
 
  bool skip_spaces (const std::string& text, size_type& pos) 
      { 
        pos = text.find_first_not_of (" \t", pos); 
        if (pos == std::string::npos) return (false); 
        return (true); 
      } 
 
  bool test_token (const std::string& text, size_type& pos, char token) 
      { 
        if (!skip_spaces (text, pos))return (false); 
 
        char c = text[pos]; 
        if (c != token) return (false); 
 
        pos++; 
        return (true); 
      } 
 
  bool parse_number (const std::string& text,  
                     size_type& pos,  
                     int& value) 
      { 
        if (pos == std::string::npos) return (false); 
         
        const char& cha = (text.at (pos)); 
        const char* ptr = &cha; 
        int items; 
 
        value = 0; 
         
        items = sscanf (ptr, "%d", &value); 
        if (items == 0)  
          { 
            return (false); 
          } 
         
        pos = text.find_first_not_of ("0123456789+- \t", pos); 
         
        return (true); 
      } 
 
  bool parse_maximum (Range::field& field, const std::string& text, size_type& pos) 
      { 
        bool result = true; 
        int maximum; 
 
        if (!skip_spaces (text, pos)) return (false); 
         
        char c = text[pos]; 
        switch (c) 
          { 
            case '0': 
            case '1': 
            case '2': 
            case '3': 
            case '4': 
            case '5': 
            case '6': 
            case '7': 
            case '8': 
            case '9': 
            case '+': 
            case '-': 
              if (!parse_number (text, pos, maximum))  
                { 
                  result = false; 
                  break; 
                } 
              else 
                { 
                  result = true; 
                  field.set_maximum (maximum); 
                } 
              break; 
            default: 
              result = false; 
              break; 
          } 
         
        return (result); 
      } 
   
  bool parse_list (Range::field& field, const std::string& text, size_type& pos) 
      { 
        bool finished = false; 
        bool result = true; 
        int value; 
 
        while (!finished) 
          { 
            if (!skip_spaces (text, pos)) return (true); 
         
            char c = text[pos]; 
            switch (c) 
              { 
              case '0': 
              case '1': 
              case '2': 
              case '3': 
              case '4': 
              case '5': 
              case '6': 
              case '7': 
              case '8': 
              case '9': 
              case '+': 
              case '-': 
                if (!parse_number (text, pos, value))  
                  { 
                    finished = true; 
                    result = false; 
                    break; 
                  } 
                else 
                  { 
                    field.add_value (value); 
                  } 
                break; 
              case ',': 
                pos++; 
                break; 
              default: 
                finished = true; 
                result = true; 
                break; 
              } 
          } 
 
        return (result); 
      } 
   
  bool parse_field (Range::field& field, const std::string& text, size_type& pos) 
      { 
        bool result = true; 
        int minimum; 
 
        if (!skip_spaces (text, pos)) return (false); 
         
        char c = text[pos]; 
        switch (c) 
          { 
            case '0': 
            case '1': 
            case '2': 
            case '3': 
            case '4': 
            case '5': 
            case '6': 
            case '7': 
            case '8': 
            case '9': 
            case '+': 
            case '-': 
              if (!parse_number (text, pos, minimum))  
                { 
                  result = false; 
                  break; 
                } 
              else 
                { 
                  field.set_minimum (minimum); 
		  result = true;
		  
                  if (test_token (text, pos, ':')) 
                    { 
			// max is optional, so we don't want to reset result here RDS 4/02
			// result = parse_maximum (field, text, pos); 

			// Optionally, look for max
			parse_maximum (field, text, pos); 
                    } 
                  else if (test_token (text, pos, ',')) 
                    { 
                      result = parse_list (field, text, pos); 
                    } 
                  else 
                    { 
                      field.set_maximum (minimum); 
                      result = true; 
                    } 
                } 
 
              break; 
            case ':': 
              pos++; 
              result = parse_maximum (field, text, pos); 
              break; 
            case '*': 
              pos++; 
              result = true; 
              break; 
            default: 
              result = false; 
              break; 
          } 
         
        return (result); 
      } 
   
  bool parse (Range& range, const std::string& text, size_type& pos) 
      { 
        bool result = true; 
        bool finished = false; 
         
        if (!skip_spaces (text, pos)) return (true); 
 
        while (!finished && pos != std::string::npos && pos < text.size()) 
          { 
            char c = text[pos]; 
         
            switch (c) 
              { 
                case '0': 
                case '1': 
                case '2': 
                case '3': 
                case '4': 
                case '5': 
                case '6': 
                case '7': 
                case '8': 
                case '9': 
                case '+': 
                case '-': 
                case ':': 
                case '*': 
                  { 
                    Range::field field; 
                    if (!parse_field (field, text, pos))  
                      { 
                        result = false; 
                        finished = true; 
                      } 
                    else 
                      { 
                        range.add (field); 
                      } 
                  } 
                  break; 
                case '/': 
                  pos++; 
                  break; 
                default: 
                  finished = true; 
                  break; 
              } 
          } 
         
        return (result); 
      } 
}; 
 
 
 
 
 
//----------------------------------------------- 
Range::field::field () 
    :
    m_minimum(0),
    m_maximum(0),
    m_indices(0),
    m_previous(0),
    m_next(0),
    m_mode(unbounded),
    m_continuation_mode(none)
{ 
} 
 
//----------------------------------------------- 
Range::field::field (const field& other) 
  : m_values (other.m_values),
    m_indexes (other.m_indexes)
{ 
  m_minimum  = other.m_minimum; 
  m_maximum  = other.m_maximum; 
  m_indices  = other.m_indices;
  m_previous = other.m_previous;
  m_next     = other.m_next;
  m_mode     = other.m_mode; 
  m_continuation_mode = other.m_continuation_mode;

} 

//----------------------------------------------- 
Range::field::field (field&& other) 
{ 
  m_minimum  = other.m_minimum; 
  m_maximum  = other.m_maximum; 
  m_values.swap (other.m_values);
  m_indexes.swap (other.m_indexes);
  m_indices  = other.m_indices;
  m_previous = other.m_previous;
  m_next     = other.m_next;
  m_mode     = other.m_mode; 
  m_continuation_mode = other.m_continuation_mode;

} 
 
//----------------------------------------------- 
Range::field::field (element_type value) 
    :
    m_minimum(value),
    m_maximum(value),
    m_indices(0),
    m_previous(0),
    m_next(0),
    m_mode(both_bounded),
    m_continuation_mode(none)
{} 
 
//----------------------------------------------- 
Range::field::field (element_type minimum, element_type maximum) 
   :
    m_indices(0),
    m_previous(0),
    m_next(0),
    m_mode(both_bounded),
    m_continuation_mode(none)
{ 
  set (minimum, maximum); 
  set_indices();    
} 
 
//----------------------------------------------- 
bool Range::field::is_valued () const  
{ 
  return (m_mode != unbounded);  
} 
 
//----------------------------------------------- 
bool Range::field::has_minimum () const 
{ 
  return ((m_mode == low_bounded) ||  
          (m_mode == both_bounded) ||  
          (m_mode == enumerated)); 
} 
 
//----------------------------------------------- 
bool Range::field::has_maximum () const 
{ 
  return ((m_mode == high_bounded) ||  
          (m_mode == both_bounded) ||  
          (m_mode == enumerated)); 
} 
 
//----------------------------------------------- 
bool Range::field::wrap_around () const 
{ 
  return (has_wrap_around == m_continuation_mode);
} 
 
 
//----------------------------------------------- 
bool 
Range::field::get_previous (element_type current, element_type& previous) const
{
  switch (m_mode) 
    { 
    case unbounded: 
      previous = current - 1; 
      if (current == std::numeric_limits<element_type>::min()) return (false);
      return (true); 
      break; 
    case low_bounded: 
      if (current == m_minimum) {
	  previous = current;
	  return (false);
      }
      previous = current - 1; 
      return (true); 
      break; 
    case high_bounded: 
      previous = current - 1; 
      if (current == std::numeric_limits<element_type>::min()) return (false);
      return (true); 
      break; 
    case both_bounded: 
      if (current == m_minimum) {
	  if (has_wrap_around == m_continuation_mode) {
	      previous = m_maximum;
	      return (true); 
	  }
	  else if (has_previous == m_continuation_mode) {
	      previous = m_previous;
	      return (true); 
	  }
	  previous = current;
	  return (false);
      }
      previous = current - 1; 
      return (true); 
      break; 
    case enumerated: 
      size_type index = get_value_index(current);
      if (index == 0) {
	  if (has_wrap_around == m_continuation_mode && m_values.size() > 0) {
	      index = m_values.size() - 1;
	      previous = m_values[index];
	      return (true); 
	  }
	  else if (has_previous == m_continuation_mode) {
	      previous = m_previous;
	      return (true); 
	  }
	  previous = current;
	  return (false);
      }
      --index;
      previous = m_values[index];
      return (true);
      break; 
    } 
 
  return (false); 
}


//----------------------------------------------- 
bool 
Range::field::get_next     (element_type current, element_type& next) const
{
  switch (m_mode) 
    { 
    case unbounded: 
      next = current + 1; 
      if (current == std::numeric_limits<element_type>::max()) return (false);
      return (true); 
      break; 
    case low_bounded: 
      next = current + 1; 
      if (current == std::numeric_limits<element_type>::max()) return (false);
      return (true); 
      break; 
    case high_bounded: 
      if (current == m_maximum) {
	  next = current;
	  return (false);
      }
      next = current + 1; 
      return (true); 
      break; 
    case both_bounded: 
      if (current == m_maximum) {
	  if (has_wrap_around == m_continuation_mode) {
	      next = m_minimum;
	      return (true); 
	  }
	  else if (has_next == m_continuation_mode) {
	      next = m_next;
	      return (true); 
	  }
	  next = current;
	  return (false);
      }
      next = current + 1; 
      return (true); 
      break; 
    case enumerated: 
      size_type index = get_value_index(current);
      if ((index == m_values.size() - 1) ||
	  (index == 0 && current != m_values[0])) {
	  if (has_wrap_around == m_continuation_mode) {
	      next = m_values[0];
	      return (true); 
	  }
	  else if (has_next == m_continuation_mode) {
	      next = m_next;
	      return (true); 
	  }
	  next = current;
	  return (false);
      }
      ++index;
      next = m_values[index];
      return (true);
      break; 
    } 
 
  return (false); 
}

//----------------------------------------------- 
bool Range::field::match_any () const 
{ 
  if (m_mode == unbounded) return (true); 
   
  return (false); 
} 
 
 
//----------------------------------------------- 
bool Range::field::overlaps_with (const field& other) const 
{ 
  typedef enum  
    { 
      done,  
      min_max,  
      max_min,  
      both, 
      both_and_enum, 
      enum_and_both, 
      all_values 
    } check_type; 
 
  static check_type check_needed[5][5] = 
  { 
    {done, done,    done,    done,          done}, 
    {done, done,    min_max, min_max,       min_max}, 
    {done, max_min, done,    max_min,       max_min}, 
    {done, max_min, min_max, both,          both_and_enum}, 
    {done, max_min, min_max, enum_and_both, all_values} 
  }; 
 
  mode other_mode = other.get_mode (); 
 
  switch (check_needed [m_mode][other_mode]) 
    { 
    case done: 
      return (true); 
    case min_max: 
      if (m_minimum <= other.get_maximum ()) return (true); 
      break; 
    case max_min: 
      if (m_maximum >= other.get_minimum ()) return (true); 
      break; 
    case both: 
      if ((m_minimum <= other.get_maximum ()) && 
          (m_maximum >= other.get_minimum ()))  
        { 
          return (true); 
        } 
      break; 
    case both_and_enum: // this is both_bounded while other is enumerated 
      if ((m_minimum <= other.get_maximum ()) && 
          (m_maximum >= other.get_minimum ()))  
        { 
          // Check if this(bb) is entirely within other(enum). 
          if ((m_minimum > other.get_minimum ()) && 
              (m_maximum < other.get_maximum ())) 
            { 
              const element_vector& ev = other.get_values (); 
              for (size_t i = 0; i < ev.size (); ++i) 
                { 
                  element_type v = ev[i]; 
                  if ((v < m_minimum) || (v > m_maximum)) return (false); 
                } 
            } 
 
          return (true); 
        } 
      break; 
    case enum_and_both: // this is enumerated while other is both_bounded 
      if ((m_minimum <= other.get_maximum ()) && 
          (m_maximum >= other.get_minimum ()))  
        { 
          // Check if other(bb) is entirely within this(enum). 
          if ((other.get_minimum () > m_minimum) && 
              (other.get_maximum () < m_maximum)) 
            { 
              const element_vector& ev = get_values (); 
              for (size_t i = 0; i < ev.size (); ++i) 
                { 
                  element_type v = ev[i]; 
                  if ((v < other.get_minimum ()) || (v > other.get_maximum ())) return (false); 
                } 
            } 
 
          return (true); 
        } 
      break; 
    case all_values: 
      // Both fields are enumerated only if there is possibility of overlap 
      if ((m_minimum <= other.get_maximum ()) && 
          (m_maximum >= other.get_minimum ()))  
        { 
          const element_vector& ev = other.get_values (); 
          for (size_t i = 0; i < m_values.size (); ++i) 
            { 
              element_type v = m_values[i]; 
              for (size_t j = 0; j < ev.size (); ++j) 
                { 
                  if (v == ev[j]) return (true); 
                } 
            } 
        } 
      break; 
    } 
 
  return (false); 
} 
 

//----------------------------------------------- 
void Range::field::clear () 
{ 
  m_minimum  = 0; 
  m_maximum  = 0; 
  m_previous = 0;
  m_next     = 0;
  m_continuation_mode = none;
  m_mode = unbounded; 
  m_values.clear (); 
} 
 
//----------------------------------------------- 
void Range::field::set (element_type minimum, element_type maximum)  
{
  m_values.clear ();

  if (minimum == maximum)
    {
      add_value (minimum);
    }
  else
    {
      m_minimum = (minimum <= maximum) ? minimum : maximum; 
      m_maximum = (maximum >= minimum) ? maximum : minimum; 
 
      m_mode = both_bounded; 
    }

  set_indices();
} 
 
//----------------------------------------------- 
void Range::field::set_minimum (element_type value)  
{ 
  if (m_mode == unbounded)  
    { 
      m_mode = low_bounded; 
      m_minimum = value; 
    } 
  else if ((m_mode == high_bounded) ||  
           (m_mode == both_bounded) ||  
           (m_mode == enumerated)) 
    { 
      set (value, get_maximum ()); 
    } 
  else 
    { 
      m_minimum = value; 
    } 

  set_indices();
} 
 
//----------------------------------------------- 
void Range::field::set_maximum (element_type value)  
{ 
  if (m_mode == unbounded)  
    { 
      m_mode = high_bounded; 
      m_maximum = value; 
    } 
  else if ((m_mode == low_bounded) ||  
           (m_mode == both_bounded) || 
           (m_mode == enumerated)) 
    { 
      set (get_minimum (), value); 
    } 
  else 
    { 
      m_maximum = value; 
    } 

  set_indices();
} 
 
//----------------------------------------------- 
void Range::field::add_value (element_type value) 
{ 
  if (m_mode == low_bounded) 
    { 
      m_values.clear (); 
      m_values.push_back (m_minimum); 
      m_maximum = m_minimum; 
      m_mode = enumerated; 
    } 
  else if (m_mode != enumerated) 
    { 
      m_values.clear (); 
      m_mode = enumerated; 
    } 
 
  for (size_type i = 0; i < m_values.size (); ++i) 
    { 
      if (m_values[i] == value) return; 
    } 
 
  m_values.push_back (value); 
  std::sort (m_values.begin (), m_values.end()); 
 
  m_minimum = m_values[0]; 
  m_maximum = m_values[m_values.size () - 1]; 

  set_indices();
} 
 
//----------------------------------------------- 
void Range::field::set (const std::vector <element_type>& values) 
{ 
  if (values.size () == 0) 
    { 
      clear (); 
      return; 
    } 
 
  for (size_type i = 0; i < values.size (); ++i) 
    { 
      add_value (values[i]); 
    } 

  set_indices();
} 
 
//----------------------------------------------- 
void Range::field::set (bool wraparound) 
{ 
    if (wraparound) {
	m_continuation_mode = has_wrap_around;
    }
} 
 
//----------------------------------------------- 
void Range::field::set_next (int next)
{
    if (has_previous == m_continuation_mode) {
	m_continuation_mode = has_both;
    }
    else {
	m_continuation_mode = has_next;
    }
    m_next = next;
}


//----------------------------------------------- 
void Range::field::set_previous (int previous)
{
    if (has_next == m_continuation_mode) {
	m_continuation_mode = has_both;
    }
    else {
	m_continuation_mode = has_previous;
    }
    m_previous = previous;
}

//----------------------------------------------- 
Range::field& Range::field::operator = (const field& other)
{
  if (this != &other) {
    m_minimum  = other.m_minimum; 
    m_maximum  = other.m_maximum; 
    m_values   = other.m_values; 
    m_indexes  = other.m_indexes;
    m_indices  = other.m_indices;
    m_previous = other.m_previous;
    m_next     = other.m_next;
    m_mode     = other.m_mode; 
    m_continuation_mode = other.m_continuation_mode;
  }

  return (*this);
}

//----------------------------------------------- 
void Range::field::operator |= (const field& other)
{
  mode other_mode = other.get_mode ();

  if (m_mode == other_mode)
    {
        /*
          x . . . .
          . x . . .
          . . x . .
          . . . x .
          . . . . x
        */
      switch (m_mode) 
        { 
          case unbounded: 
            break; 
          case high_bounded: 
            if (other.get_maximum () > m_maximum) m_maximum = other.get_maximum ();
            break; 
          case low_bounded: 
            if (other.get_minimum () < m_minimum) m_minimum = other.get_minimum ();
            break; 
          case enumerated:
          {
            const element_vector& ev = other.get_values ();

            for (size_t i = 0; i < ev.size (); ++i) 
              { 
                add_value (ev[i]);
              } 
          }
            break; 
          default:  // both_bounded 
              /**
               *  If there is no overlap we should build a multi-segment specification.
               *  The current algorithm is only correct if the overlap in not empty !!
               *   A multi-segment specification might also be implemented as an 
               *  expanded enumerated set (not very optimized !!)
               */
            if (other.get_maximum () > m_maximum) m_maximum = other.get_maximum ();
            if (other.get_minimum () < m_minimum) m_minimum = other.get_minimum ();

            break; 
        } 
    }
  else if ((m_mode == unbounded) || (other_mode == unbounded))
    {
        /*
          o x x x x
          x o . . .
          x . o . .
          x . . o .
          x . . . o
         */
      clear ();
    }
  else if ((m_mode == low_bounded) && (other_mode == high_bounded))
    {
        /**
         *  If there is no overlap we should build a multi-segment specification.
         *  The current algorithm is only correct if the overlap in not empty !!
         *
         *   (in addition, the expanded solution - to enumerated - is not possible
         *    due to the unbounded nature of this mode)
         */


        /*
          o o o o o
          o o x . .
          o . o . .
          o . . o .
          o . . . o
        */
      clear ();
    }
  else if ((m_mode == high_bounded) && (other_mode == low_bounded))
    {
        /**
         *  If there is no overlap we should build a multi-segment specification.
         *  The current algorithm is only correct if the overlap in not empty !!
         *
         *   (in addition, the expanded solution - to enumerated - is not possible
         *    due to the unbounded nature of this mode)
         */


        /*
          o o o o o
          o o o . .
          o x o . .
          o . . o .
          o . . . o
         */
      clear ();
    }
  else
    {
        // all other cases...

      if (has_minimum () && other.has_minimum ())
        {
            /*
              o o o o o
              o o o x x
              o o o . .
              o x . o x
              o x . x o
            */

          if (other.get_minimum () < m_minimum) set_minimum (other.get_minimum ());
        }

      if (has_maximum () && other.has_maximum ())
        {
            /*
              o o o o o
              o o o . .
              o o o x x
              o . x o x
              o . x x o
            */

          if (other.get_maximum () > m_maximum) set_maximum (other.get_maximum ());
        }
    }

  set_indices();
}

//----------------------------------------------- 
Range::field::operator std::string () const 
{ 
  std::string result; 
  char temp[20]; 

  if (!is_valued ()) 
    { 
      result += "*"; 
    } 
  else  
    { 
      element_type minimum = get_minimum (); 
      element_type maximum = get_maximum (); 
 
      if (!has_maximum ()) 
        { 
          sprintf (temp, "%d", minimum); 
          result += temp; 
          result += ":"; 
        } 
      else if (!has_minimum ()) 
        { 
          sprintf (temp, "%d", maximum); 
          result += ":"; 
          result += temp; 
        } 
      else if (minimum == maximum) 
        { 
          sprintf (temp, "%d", minimum); 
          result += temp; 
        } 
      else 
        { 
          if (get_mode () == field::enumerated) 
            { 
              for (size_type i = 0; i < get_indices (); ++i) 
                { 
                  if (i > 0) result += ","; 
                  sprintf (temp, "%d", get_value_at (i)); 
                  result += temp; 
                } 
            } 
          else 
            { 
              sprintf (temp, "%d", minimum); 
              result += temp; 
              sprintf (temp, "%d", maximum); 
              result += ":"; 
                  result += temp; 
            } 
        } 
    } 
 
  return (result); 
} 

//-----------------------------------------------
bool 
Range::field::operator == (const field& other) const
{
    if (m_mode != other.m_mode) 	return false;
    if (m_minimum != other.m_minimum) 	return false;
    if (m_maximum != other.m_maximum) 	return false;
    if (m_values  != other.m_values) 	return false;
    return (true);
}

//----------------------------------------------- 
bool 
Range::field::operator != (const field& other) const
{
    return (!((*this) == other));
}
 

//----------------------------------------------- 
void Range::field::show() const
{
    
  std::cout << "min/max " << m_minimum << " " << m_maximum << " "; 
  std::cout << "values  ";
  for (size_type i = 0; i < m_values.size(); ++i) {
      std::cout << m_values[i] << " ";
  }
  std::cout << "indexes  ";
  for (size_type i = 0; i < m_indexes.size(); ++i) {
      std::cout << m_indexes[i] << " ";
  }
  std::cout << "indices  " << m_indices << " ";
  std::cout << "prev  " << m_previous << " ";
  std::cout << "next  " << m_next << " ";

  std::cout << "mode  ";
  switch (m_mode) { 
  case Range::field::unbounded: 
      std::cout << "unbounded  ";
      break; 
  case Range::field::low_bounded: 
      std::cout << "low_bounded  ";
      break; 
  case Range::field::high_bounded: 
      std::cout << "high_bounded  ";
      break; 
  case Range::field::both_bounded: 
      std::cout << "both_bounded  ";
      break; 
  case Range::field::enumerated: 
      std::cout << "enumerated  ";
      break; 
  } 

  std::cout << "cont mode  ";
  switch (m_continuation_mode) { 
  case Range::field::none: 
      std::cout << "none  ";
      break; 
  case Range::field::has_next: 
      std::cout << "has_next  ";
      break; 
  case Range::field::has_previous: 
      std::cout << "has_previous  ";
      break; 
  case Range::field::has_both:
      std::cout << "has_both  ";
      break; 
  case Range::field::has_wrap_around:
      std::cout << "has_wrap_around  ";
      break; 
  }
  std::cout << std::endl;
}



//----------------------------------------------- 
void Range::field::optimize()
{

    /// Check mode - switch from enumerated to both_bounded if possible
    check_for_both_bounded();

    /// Create index table from value table
    create_index_table();

}

//----------------------------------------------- 
void Range::field::set_indices()
{
    /// Set the number of indices
    m_indices = 1;
    if (m_mode == both_bounded) { 
	m_indices = m_maximum - m_minimum + 1; 
    } 
    else if (m_mode == enumerated) { 
	m_indices = m_values.size (); 
    } 
} 

//----------------------------------------------- 
void Range::field::check_for_both_bounded()
{
    if (m_mode == enumerated && m_values.size() > 0) {
	element_type last = m_values[0];	
	for (size_type i = 1; i < m_values.size (); ++i) { 
	    if (m_values[i] > last + 1) return;
	    last = m_values[i];
        }

	// Is both bounded - swith mode
	m_minimum = m_values[0];
	m_maximum = m_values[m_values.size() - 1];

//  	for (size_type i = 0; i < m_values.size (); ++i) { 
//  	    std::cout << m_values[i] << " ";
//          }
//  	std::cout << " min/max " << m_minimum << " " << m_maximum << std::endl;

	m_mode = both_bounded; 
	m_values.clear (); 
    }
}

//----------------------------------------------- 
void Range::field::create_index_table()
{
    /// Create index table from value table
    if (m_mode == enumerated && m_values.size() > 0) {
	size_type size = m_maximum - m_minimum + 1;
	// return if we are over the maximum desired vector table size	
	if (size > max_indexes) {
	    m_indexes.clear();
	    return;
	}
	// Set up vectors for decoding
	m_indexes = std::vector<size_type>(size, 0);
	size_type index = 0;
	for (size_type i = 0; i < m_values.size(); ++i) {
	    if ((m_values[i]- m_minimum) < (int)size) {
		m_indexes[m_values[i] - m_minimum] = index;
		index++;
	    }
	    else {
		std::cout << "size, value, index, i " 
			  << size << " " << m_values[i] << " "
			  << index << " " << i << "  min, max " 
			  << m_minimum << " " 
			  << m_maximum 
			  << std::endl;
	    }
	}
    }
}

 
    // Constructors 
//----------------------------------------------- 
Range::Range () 
{ 
} 
 
//----------------------------------------------- 
Range::Range (const Range& other) 
{ 
  m_fields = other.m_fields; 
} 
 
//----------------------------------------------- 
Range& Range::operator= (const Range& other) 
{ 
  if (this != &other)
    m_fields = other.m_fields; 
  return *this;
} 
 
//----------------------------------------------- 
Range::Range (Range& other, bool) 
{ 
  m_fields.swap (other.m_fields);
} 
 
//----------------------------------------------- 
Range::Range (const Range& other, size_type start) 
{ 
  if (start < other.fields ()) 
    { 
      field_vector::const_iterator it = other.m_fields.begin (); 
      it += start; 
 
      m_fields.insert (m_fields.end (), it, other.m_fields.end ()); 
    } 
} 
 
//----------------------------------------------- 
//  Range::Range (const std::string& text) 
//  { 
//    build (text); 
//  } 
 
//----------------------------------------------- 
Range::Range (const ExpandedIdentifier& root) 
{ 
  // Construct from a root (i.e. add wild card for below) 
  build (root); 
} 
 
//----------------------------------------------- 
void Range::build (const std::string& text) 
{ 
  RangeParser parser; 
 
  parser.run (*this, text); 
} 
 
//----------------------------------------------- 
void Range::build (const ExpandedIdentifier& root) 
{ 
    // Construct from a root  
  m_fields.clear (); 
 
  for (size_type i = 0; i < root.fields (); ++i) 
    { 
      m_fields.push_back (field (root[i])); 
    } 
} 
 
    // Modifications 
//----------------------------------------------- 
void Range::add () 
{ 
  m_fields.push_back (field ()); 
} 
 
//----------------------------------------------- 
void Range::add (element_type value) 
{ 
  m_fields.push_back (field (value)); 
} 
 
//----------------------------------------------- 
void Range::add (element_type minimum, element_type maximum) 
{ 
  m_fields.push_back (field (minimum, maximum)); 
} 
 
//----------------------------------------------- 
void Range::add_minimum (element_type minimum) 
{ 
  field f; 
   
  f.set_minimum (minimum); 
  m_fields.push_back (f); 
} 
 
//----------------------------------------------- 
void Range::add_maximum (element_type maximum) 
{ 
  field f; 
 
  f.set_maximum (maximum); 
  m_fields.push_back (f); 
} 
 
/// Add a range specified using a field  
void Range::add (const field& f) 
{
  m_fields.emplace_back (f); 
} 
 
/// Add a range specified using a field, using move semantics.
void Range::add (field&& f) 
{
  m_fields.emplace_back (std::move(f));
} 
 
/// Append a subrange 
void Range::add (const Range& subrange) 
{ 
  for (size_t i = 0; i < subrange.fields (); ++i) 
    { 
      const field& f = subrange[i]; 
      m_fields.push_back (f); 
    } 
} 

void Range::add (Range&& subrange) 
{
  if (m_fields.empty())
    m_fields.swap (subrange.m_fields);
  else {
    size_t sz = subrange.m_fields.size();
    m_fields.reserve (m_fields.size() + sz);
    for (size_t i = 0; i < sz; ++i) 
    { 
      m_fields.emplace_back (std::move(subrange.m_fields[i]));
    }
  }
} 

//----------------------------------------------- 
void Range::clear () 
{ 
  m_fields.clear (); 
} 
 
//----------------------------------------------- 
int Range::match (const ExpandedIdentifier& id) const 
{ 
  size_type my_fields = m_fields.size (); 
  const size_type id_fields = id.fields (); 
 
    // Remove trailing wild cards since they are meaningless. 
  while ((my_fields > 1) && 
         (!m_fields[my_fields-1].is_valued ())) 
    { 
      my_fields--; 
    } 
 
    // Ranges with only wild cards always match. 
  if (my_fields == 0) return (1); 
 
    // More fields in the range than in the identifier will never match. 
  //if (my_fields > id_fields) return (0); 

  // Allow match for id shorter than range - assume "wildcards" for
  // missing id fields
  size_type nfields =  (my_fields > id_fields) ? id_fields : my_fields;
 
    // Test fields one by one. 
  //for (size_type field_number = 0; field_number < my_fields; field_number++) 
  for (size_type field_number = 0; field_number < nfields; field_number++) 
    { 
      const field& f = m_fields[field_number]; 
 
      if (!f.match (id[field_number])) return (0); 
    } 
 
    // All conditions match. 
  return (1); 
} 
 
// Accessors 
 
//----------------------------------------------- 
ExpandedIdentifier Range::minimum () const 
{ 
  size_type my_fields = m_fields.size (); 
  ExpandedIdentifier result; 
   
    // Remove trailing wild cards since they are meaningless. 
  while ((my_fields > 1) && 
         (!m_fields[my_fields-1].has_minimum ()))  
    { 
      my_fields--; 
    } 
   
    // Ranges with only wild cards: set first field of min to 0 
  if (my_fields == 0) {
    result << 0;
    return result; // Don't combine these two lines --- it inhibits RVO.
  }
   
    // Copy fields to result - look for wild cards 
  for (size_type field_number = 0; field_number < my_fields; field_number++)  
    { 
      const field& f = m_fields[field_number]; 
       
      if (!f.has_minimum ())  
        { 
            // Wilds card -> set field to 0 
          result << 0; 
        }        
      else  
        { 
            // Valued field 
          result << f.get_minimum (); 
        } 
    } 
   
  return (result); 
} 
 
//----------------------------------------------- 
ExpandedIdentifier Range::maximum () const 
{ 
  size_type my_fields = m_fields.size (); 
  ExpandedIdentifier result; 
   
    // Remove all by the last trailing wild card, extra ones are 
    // meaningless. 
  while ((my_fields > 1) && 
         (!m_fields[my_fields-1].has_maximum ())) 
    { 
      my_fields--; 
    } 
 
    // Ranges with only wild cards: set first field of min to ExpandedIdentifier::max_value 
  if (my_fields == 0) {
    result << ExpandedIdentifier::max_value;
    return result; // Don't combine these two lines --- it inhibits RVO.
  }
 
    // Copy fields to result - look for wild cards 
  for (size_type field_number = 0; field_number < my_fields; field_number++)  
    { 
      const field& f = m_fields[field_number]; 
 
      if (!f.has_maximum ())  
        { 
            // Wilds card  
          if (field_number == 0) 
            { 
                // For 1st field set it to ExpandedIdentifier::max_value 
              result << ExpandedIdentifier::max_value; 
            } 
          else 
            { 
                // 
                // For subsequent fields, set do ++ for field-1 
                // This requires rebuilding the result 
                // 
              ExpandedIdentifier new_result; 
 
              for (size_type new_field_number = 0;  
                   new_field_number < (field_number - 1); 
                   ++new_field_number) 
                { 
                  new_result << result[new_field_number]; 
                } 
 
              element_type last = result[field_number - 1]; 
 
              new_result << ((ExpandedIdentifier::max_value == last) ? last : last + 1); 
              new_result << 0; 
 
              assert ( result.fields () == new_result.fields () ); 
 
              result = new_result; 
            } 
        } 
      else 
        { 
            // Normal field 
          result << f.get_maximum (); 
        } 
    } 
   
    return (result); 
} 
 
Range::size_type Range::cardinality () const 
{ 
  size_type result = 1; 
 
  const Range& me = *this; 
 
  for (size_type i = 0; i < fields (); ++i) 
    { 
      const field& f = me[i]; 
 
      result *= f.get_indices (); 
    } 
 
  return (result); 
} 
 

/**
 *   Get the cardinality from the beginning up to the given ExpandedIdentifier
 */
Range::size_type Range::cardinalityUpTo (const ExpandedIdentifier& id) const 
{ 
  size_type result = 0; 

//    std::cout << " Range::cardinality: id, fields " << (std::string) id 
//  	    << " " << id.fields() << " " << fields() << " "
//  	    << (std::string)(*this)  << std::endl;

  if (id.fields() != fields()) return (result);
 
  const Range& me = *this; 

  // Check if we are above or below this range
  if (id < minimum()) return 0;
  if (maximum() < id) return cardinality();
  
  // Collect the indices of a match
  std::vector<size_type> indices(id.fields (), 0);
  bool is_match = true;
  size_type level = 0;
  for (; level < id.fields (); ++level) {

      const field& f = me[level]; 

      // Require all fields to be bounded or enumerated
      if (!(f.get_mode() == Range::field::both_bounded ||
	    f.get_mode() == Range::field::enumerated)) return 0;
      
      //int index = 0;  // Contains number of nodes below match

      if (f.get_mode() == Range::field::enumerated) {
	  // Continue testing for a match
	  size_type max = f.get_values().size() - 1;
	  if (f.get_values()[max] < id[level]) {
	      // above max
	      is_match = false;
	      indices[level] = max + 1;
	  }
	  else {
	      for (size_type j = 0; j < f.get_values().size(); ++j) {
		  if (id[level] <= f.get_values()[j]) {
		      if (id[level] != f.get_values()[j]) {
			  // between two values or below first one
			  is_match = false;
		      }
		      indices[level] = j;
		      break;
		  }
	      }
	  }
      }
      else {
	  if (f.get_maximum() < id[level]) {
	      // above max
	      is_match = false;
	      indices[level] = f.get_maximum() - f.get_minimum() + 1;
	  }
	  else if (id[level] < f.get_minimum()) {
	      // below min
	      is_match = false;
	      indices[level] = 0;
	  }
	  else {
	      indices[level] = id[level] - f.get_minimum();
	  }
      }

//        std::cout << " level, id, field, indices[[level], is enum, is_match " << level << " " 
//  		<< id[level] << " " << (std::string)f << " " 
//  		<< indices[level] << " " << (f.get_mode() == Range::field::enumerated) << " " 
//  		<< is_match << std::endl;


      if (!is_match) break;
  }
  
  // Calculate the cardinality
  if (level < id.fields ()) ++level;
  for (size_type j = 0; j < level; ++j) {
      size_type card = indices[j];
      for (size_type k = j + 1; k < id.fields(); ++k) {

	  const field& f = me[k]; 

	  card *= f.get_indices();
      }
      result += card;

//        std::cout << " j, indices, card " << j << " " 
//  		<< indices[j] << " " << card << std::endl;
  }
  
      
  return result;
  
} 
 
 
/** 
 *   Check overlap between two Ranges : 
 * 
 *   As soon as one pair of corresponding fields do not match, 
 *   the global overlap is empty. 
 */ 
bool Range::overlaps_with (const Range& other) const 
{ 
  const Range& me = *this; 
 
  if ((fields () == 0) || (other.fields () == 0)) return (false); 
 
  for (size_type i = 0; i < std::min (fields (), other.fields ()); ++i) 
    { 
      const field& f1 = me[i]; 
      const field& f2 = other[i]; 
 
      if (!f1.overlaps_with (f2)) return (false); 
    } 
 
  return (true); 
} 
 
//----------------------------------------------- 
void Range::show () const 
{
  show (std::cout);
}

void Range::show (std::ostream& s) const 
{ 
  const Range& me = *this; 
 
  s << (std::string) me << " (";
 
  int allbits = 0;

  for (size_type i = 0; i < fields (); ++i) 
    { 
      const field& f = me[i]; 

      if (i > 0) s << "+";

      size_type indices = f.get_indices ();

      int bits = 1;
      indices--;
      if (indices > 0)
        {
          bits = 0;
          while (indices > 0)
            {
              indices /= 2;
              bits++;
            }
        }

      allbits += bits;
 
      s << bits; 
    } 

  s << "=" << allbits << ") ";
} 
 
//----------------------------------------------- 
Range::operator std::string () const 
{ 
  std::string result; 
 
  size_type my_fields = m_fields.size (); 

    // Remove trailing wild cards since they are meaningless. 
  while ((my_fields > 1) && 
         (!m_fields[my_fields-1].is_valued ())) 
    { 
      my_fields--; 
    } 
 
  if (my_fields == 0) return (result); 
 
    // print fields one by one. 
  for (size_type field_number = 0; field_number < my_fields; field_number++) 
    { 
      const field& f = m_fields[field_number]; 

      if (field_number > 0) result += "/"; 
      result += (std::string) f; 
    } 
 
  return (result); 
} 
 
//----------------------------------------------- 
bool 
Range::operator == (const Range& other) const
{
    if (m_fields.size() != other.m_fields.size()) return false;
    field_vector::const_iterator it1  = m_fields.begin();
    field_vector::const_iterator it2  = other.m_fields.begin();
    field_vector::const_iterator last = m_fields.end();
    for (; it1 != last; ++it1, ++it2) {
	if ((*it1) != (*it2)) return false;
    }
    return (true);
}


//----------------------------------------------- 
bool 
Range::operator != (const Range& other) const
{
    return (!((*this) == other));
}



//----------------------------------------------- 
Range::identifier_factory Range::factory_begin () 
{ 
  const Range& me = *this; 
 
  return (identifier_factory (me)); 
} 
 
//----------------------------------------------- 
Range::const_identifier_factory Range::factory_begin () const 
{ 
  const Range& me = *this; 
 
  return (const_identifier_factory (me)); 
} 
 
//----------------------------------------------- 
Range::identifier_factory Range::factory_end () 
{ 
  static Range r; 
  static identifier_factory factory (r); 
 
  return (factory); 
} 
 
//----------------------------------------------- 
Range::const_identifier_factory Range::factory_end () const 
{ 
  static const_identifier_factory factory; 
 
  return (factory); 
} 
 
//----------------------------------------------- 
Range::identifier_factory::identifier_factory () : m_range (0) 
{ 
} 
 
//----------------------------------------------- 
Range::identifier_factory::identifier_factory (const Range& range) : m_range (&range) 
{ 
  /** 
   *    Fill all running identifiers 
   *    m_id : the current id 
   *    m_min : the set of low bounds 
   *    m_max : the set of high bounds 
   */ 
  for (Range::size_type i = 0; i < range.fields (); ++i) 
    { 
      element_type minimum; 
      element_type maximum; 
 
      m_indices.push_back (0); 
 
      const field& f = range[i]; 
 
      switch (f.get_mode ()) 
        { 
        case Range::field::unbounded: 
          m_id << 0; 
          m_min << 0; 
          m_max << 0; 
          break; 
        case Range::field::low_bounded: 
          minimum = f.get_minimum (); 
          m_id << minimum; 
          m_min << minimum; 
          m_max << minimum; 
          break; 
        case Range::field::high_bounded: 
          maximum = f.get_maximum (); 
          m_id << maximum; 
          m_min << maximum; 
          m_max << maximum; 
          break; 
        case Range::field::both_bounded: 
        case Range::field::enumerated: 
          minimum = f.get_minimum (); 
          maximum = f.get_maximum (); 
          m_id << minimum; 
          m_min << minimum; 
          m_max << maximum; 
          break; 
        } 
    } 
} 
 
//----------------------------------------------- 
Range::identifier_factory& Range::identifier_factory::operator = (const identifier_factory& other) 
{ 
  if (&other != this) {
    m_indices = other.m_indices; 
    m_id      = other.m_id; 
    m_min     = other.m_min; 
    m_max     = other.m_max; 
    m_range   = other.m_range; 
  }
 
  return (*this); 
} 
 
//----------------------------------------------- 
void Range::identifier_factory::operator ++ () 
{ 
  if (m_id.fields () == 0) return; 
 
  size_type fields = m_id.fields (); 
  size_type i = fields - 1; 
 
    // 
    // Starting from the end, we try to increment the m_id fields 
    // If at a given position it's not possible (max reached) 
    // then we move back one pos and try again. 
    // 
    //  As soon as increment is possible, then the rest of the m_id 
    // is reset to min values. 
    // 
 
  for (;;) 
    { 
      const field& f = (*m_range)[i]; 
      bool done = false; 
      element_type value = 0; 
 
      if (f.get_mode () == Range::field::enumerated) 
        { 
          Range::size_type index = m_indices[i]; 
          index++; 
          if (index < f.get_indices ()) 
            { 
              m_indices[i] = index; 
              value = f.get_value_at (index); 
              done = true; 
            } 
        } 
      else 
        { 
          value = m_id[i]; 
 
          if (value < m_max[i]) 
            { 
              /** 
               *   The local range is not exceeded. 
               *   increase the value then reset the remaining fields. 
               */ 
              ++value; 
              done = true; 
            } 
        } 
 
      if (done) 
        { 
          m_id[i] = value; 
           
          for (++i; i < fields; ++i) 
            { 
              m_indices[i] = 0; 
              m_id[i] = m_min[i]; 
            } 
           
          break; 
        } 
 
      /** 
       *  The current range field was exhausted 
       *  check the previous one. 
       */ 
 
      if (i == 0)  
        { 
          m_id.clear (); 
          break; 
        } 
      else 
        { 
          --i; 
        } 
    } 
} 
 
//----------------------------------------------- 
const ExpandedIdentifier& Range::identifier_factory::operator * () const 
{ 
  return (m_id); 
} 
 
//----------------------------------------------- 
bool Range::identifier_factory::operator == (const identifier_factory& other) const 
{ 
  if (m_id == other.m_id) return (true); 
  return (false); 
} 
 
//----------------------------------------------- 
bool Range::identifier_factory::operator != (const identifier_factory& other) const 
{ 
  if (m_id != other.m_id) return (true); 
  return (false); 
} 
 
//----------------------------------------------- 
Range::const_identifier_factory::const_identifier_factory () : m_range (0) 
{ 
} 
 
//----------------------------------------------- 
Range::const_identifier_factory::const_identifier_factory (const Range& range) :  
  m_range (&range) 
{ 
  /** 
   *    Fill all running identifiers 
   *    m_id : the current id 
   *    m_min : the set of low bounds 
   *    m_max : the set of high bounds 
   */ 
  for (Range::size_type i = 0; i < range.fields (); ++i) 
    { 
      element_type minimum; 
      element_type maximum; 
 
      m_indices.push_back (0); 
 
      const field& f = range[i]; 
 
      switch (f.get_mode ()) 
        { 
        case Range::field::unbounded: 
          m_id << 0; 
          m_min << 0; 
          m_max << 0; 
          break; 
        case Range::field::low_bounded: 
          minimum = f.get_minimum (); 
          m_id << minimum; 
          m_min << minimum; 
          m_max << minimum; 
          break; 
        case Range::field::high_bounded: 
          maximum = f.get_maximum (); 
          m_id << maximum; 
          m_min << maximum; 
          m_max << maximum; 
          break; 
        case Range::field::both_bounded: 
        case Range::field::enumerated: 
          minimum = f.get_minimum (); 
          maximum = f.get_maximum (); 
          m_id << minimum; 
          m_min << minimum; 
          m_max << maximum; 
          break; 
        } 
    } 
} 
 
//----------------------------------------------- 
Range::const_identifier_factory& Range::const_identifier_factory::operator = (const const_identifier_factory& other) 
{ 
  if (&other != this) {
    m_indices = other.m_indices; 
    m_id      = other.m_id; 
    m_min     = other.m_min; 
    m_max     = other.m_max; 
    m_range   = other.m_range; 
  }
 
  return (*this); 
} 
 
//----------------------------------------------- 
void Range::const_identifier_factory::operator ++ () 
{ 
  if (m_id.fields () == 0) return; 
 
  size_type fields = m_id.fields (); 
  size_type i = fields; 
 
    // 
    // Starting from the end, we try to increment the m_id fields 
    // If at a given position it's not possible (max reached) 
    // then we move back one pos and try again. 
    // 
    //  As soon as increment is possible, then the rest of the m_id 
    // is reset to min values. 
    // 
 
  for (;;) 
    { 
      const field& f = (*m_range)[i]; 
      bool done = false; 
      element_type value = 0; 
 
      if (f.get_mode () == Range::field::enumerated) 
        { 
          Range::size_type index = m_indices[i]; 
          index++; 
          if (index < f.get_indices ()) 
            { 
              m_indices[i] = index; 
              value = f.get_value_at (index); 
              done = true; 
            } 
        } 
      else 
        { 
          value = m_id[i]; 
 
          if (value < m_max[i]) 
            { 
              /** 
               *   The local range is not exceeded. 
               *   increase the value then reset the remaining fields. 
               */ 
              ++value; 
              done = true; 
            } 
        } 
 
      if (done) 
        { 
          m_id[i] = value; 
           
          for (++i; i < fields; ++i) 
            { 
              m_indices[i] = 0; 
              m_id[i] = m_min[i]; 
            } 
           
          break; 
        } 
 
      /** 
       *  The current range field was exhausted 
       *  check the previous one. 
       */ 
 
      if (i == 0)  
        { 
          m_id.clear (); 
          break; 
        } 
      else 
        { 
          --i; 
        } 
    } 
} 
 
//----------------------------------------------- 
const ExpandedIdentifier& Range::const_identifier_factory::operator * () const 
{ 
  return (m_id); 
} 
 
//----------------------------------------------- 
bool Range::const_identifier_factory::operator == (const const_identifier_factory& other) const 
{ 
  if (m_id == other.m_id) return (true); 
  return (false); 
} 
 
//----------------------------------------------- 
bool Range::const_identifier_factory::operator != (const const_identifier_factory& other) const 
{ 
  if (m_id != other.m_id) return (true); 
  return (false); 
} 
 
 
 
class MultiRangeParser 
{ 
public: 
  typedef MultiRange::size_type size_type; 
 
  MultiRangeParser () 
      { 
        m_multirange = 0; 
      } 
 
  bool run (const std::string& text, MultiRange& multirange) 
      { 
        m_multirange = &multirange; 
        multirange.clear (); 
        size_type pos = 0; 
        return (parse (text, pos)); 
      } 
 
private: 
 
  bool parse_number (const std::string& text,  
                     size_type& pos,  
                     int& value) 
      { 
        if (pos == std::string::npos) return (false); 
         
        const char& cha = (text.at (pos)); 
        const char* ptr = &cha; 
        int items; 
 
        value = 0; 
         
        items = sscanf (ptr, "%d", &value); 
        if (items == 0)  
          { 
            return (false); 
          } 
         
        pos = text.find_first_not_of ("0123456789+- \t", pos); 
         
        return (true); 
      } 
 
  bool parse_field (const std::string& text, size_type& pos) 
      { 
        bool result = true; 
        int minimum; 
        int maximum; 
 
        if (!skip_spaces (text, pos)) return (false); 
         
        char c = text[pos]; 
        switch (c) 
          { 
            case '0': 
            case '1': 
            case '2': 
            case '3': 
            case '4': 
            case '5': 
            case '6': 
            case '7': 
            case '8': 
            case '9': 
            case '+': 
            case '-': 
              if (true) 
                { 
                  // 
                  // The current range is considered 
                  // 
                  Range& r = m_multirange->back (); 
                   
                  if (!parse_number (text, pos, minimum))  
                    { 
                      result = false; 
                      break; 
                    } 
 
                  if (test_token (text, pos, ':')) 
                    { 
                      if (!parse_number (text, pos, maximum))  
                        { 
                          r.add_minimum ((MultiRange::element_type) minimum); 
                        } 
                      else 
                        { 
                          r.add ((MultiRange::element_type) minimum,  
                                 (MultiRange::element_type) maximum); 
                        } 
                    } 
                  else 
                    { 
                      r.add ((MultiRange::element_type) minimum); 
                    } 
                } 
 
              break; 
            case ':': 
              pos++; 
              if (true) 
                { 
                  Range& r = m_multirange->back (); 
 
                  if (!parse_number (text, pos, maximum))  
                    { 
                      result = false; 
                    } 
                  else 
                    { 
                      r.add_maximum ((MultiRange::element_type) maximum); 
                    } 
                } 
 
              break; 
            case '*': 
              pos++; 
              if (true) 
                { 
                  Range& r = m_multirange->back (); 
                  r.add (); 
                } 
 
              break; 
            default: 
              result = false; 
          } 
         
        return (result); 
      } 
   
  bool parse_fields (const std::string& text, size_type& pos) 
      { 
        bool finished = false; 
        bool result = true; 
         
        while (!finished) 
          { 
            if (!skip_spaces (text, pos)) 
              { 
                result = false; 
                break; 
              } 
             
            char c = text[pos]; 
             
            switch (c) 
              { 
                case '0': 
                case '1': 
                case '2': 
                case '3': 
                case '4': 
                case '5': 
                case '6': 
                case '7': 
                case '8': 
                case '9': 
                case '+': 
                case '-': 
                case ':': 
                case '*': 
                  if (!parse_field (text, pos))  
                    { 
                      result = false; 
                      finished = true; 
                    } 
                  break; 
                case '/': 
                  pos++; 
                  break; 
                default: 
                  finished = true; 
                  break; 
              } 
          } 
         
        return (result); 
      } 
   
  bool parse (const std::string& text, size_type& pos) 
      { 
        bool result = true; 
        bool finished = false; 
         
        if (!skip_spaces (text, pos)) return (true); 
 
        while (!finished) 
          { 
            char c = text[pos]; 
         
            switch (c) 
              { 
                case '0': 
                case '1': 
                case '2': 
                case '3': 
                case '4': 
                case '5': 
                case '6': 
                case '7': 
                case '8': 
                case '9': 
                case '+': 
                case '-': 
                case ':': 
                case '*': 
                  m_multirange->add_range (); 
                  if (!parse_fields (text, pos))  
                    { 
                      result = false; 
                      finished = true; 
                    } 
                  break; 
                case '|': 
                  pos++; 
                  break; 
                default: 
                  finished = true; 
                  break; 
              } 
          } 
         
        return (result); 
      } 
   
  bool skip_spaces (const std::string& text, size_type& pos) 
      { 
        pos = text.find_first_not_of (" \t", pos); 
        if (pos == std::string::npos) return (false); 
        return (true); 
      } 
 
  bool test_token (const std::string& text, size_type& pos, char token) 
      { 
        if (!skip_spaces (text, pos))return (false); 
 
        char c = text[pos]; 
        if (c != token) return (false); 
 
        pos++; 
        return (true); 
      } 
 
  MultiRange* m_multirange; 
}; 
 
 
//----------------------------------------------- 
MultiRange::MultiRange () 
    :
    m_it_count(0)
{ 
} 
 
//----------------------------------------------- 
MultiRange::MultiRange (const MultiRange& other) 
    : 
    m_ranges (other.m_ranges),
    m_it_count(0)
{ 
} 

//----------------------------------------------- 
MultiRange& MultiRange::operator= (const MultiRange& other) 
{
  if (this != &other) {
    if (m_it_count > 0) {
      // Can't do this if we still have iterators referencing us.
      std::abort();
    }
    m_ranges = other.m_ranges;
    m_it_count = other.m_it_count;
  }
  return *this;
}

/** 
 *    The reduction algorithm is no longer provided, since it does not 
 *   work in the most general case. 
 *    See previous and temptative implementation in the RangeReduction.cxx  
 *   source file. 
 */ 
MultiRange::MultiRange (const Range& r, const Range& s) 
    : 
    m_it_count(0)
{ 
  m_ranges.push_back (r); 
  m_ranges.push_back (s); 
} 
 
//----------------------------------------------- 
MultiRange::MultiRange (const std::string& text) 
    : 
    m_it_count(0)
{ 
  build (text); 
} 
 
//----------------------------------------------- 
void MultiRange::build (const std::string& text) 
{ 
  static MultiRangeParser parser; 
  parser.run (text, *this); 
} 
 
//----------------------------------------------- 
void MultiRange::clear () 
{ 
  m_ranges.clear (); 
} 
 
//----------------------------------------------- 
void MultiRange::add (const Range& range) 
{ 
    // Add new range ONLY if an equivalent does NOT exist
    for (size_type i = 0; i < m_ranges.size(); ++i) {
	const Range& test_range = m_ranges[i];
	if (test_range == range) return;
    }
    m_ranges.push_back (range); 
} 
 
//----------------------------------------------- 
void MultiRange::add (Range& range, bool) 
{ 
    // Add new range ONLY if an equivalent does NOT exist
    for (size_type i = 0; i < m_ranges.size(); ++i) {
	const Range& test_range = m_ranges[i];
	if (test_range == range) return;
    }
    m_ranges.emplace_back (range, true);
} 
 
//----------------------------------------------- 
void MultiRange::add (const ExpandedIdentifier& id) 
{ 
  m_ranges.push_back (Range (id)); 
} 
 
//----------------------------------------------- 
void MultiRange::remove_range (const ExpandedIdentifier& id) 
{ 
    // Remove all ranges for which id matches
    range_vector::iterator first=m_ranges.begin();
    range_vector::iterator last=m_ranges.end();
    while (first != last) {
	if (first->match(id)) {
	    m_ranges.erase(first);  // remove range
	    first=m_ranges.begin(); // restart loop every time
	    last=m_ranges.end();    // when range removed 
	} else {
	    ++first;
	}
    }
}

 
//----------------------------------------------- 
Range& MultiRange::add_range () 
{ 
  size_type size = m_ranges.size (); 
  m_ranges.resize (size + 1); 
  return (m_ranges.back ()); 
} 
 
//----------------------------------------------- 
Range& MultiRange::back () 
{ 
  return (m_ranges.back ()); 
} 
 
//----------------------------------------------- 
int MultiRange::match (const ExpandedIdentifier& id) const 
{ 
  range_vector::size_type i; 
 
  for (i = 0; i < m_ranges.size (); ++i) 
    { 
      const Range& r = m_ranges[i]; 
 
      if (r.match (id)) return (1); 
    } 
 
  return (0); 
} 
 
//----------------------------------------------- 
const Range& MultiRange::operator [] (MultiRange::size_type index) const 
{ 
  static const Range null_range; 
 
  if (index >= m_ranges.size ()) return (null_range); 
 
  return (m_ranges[index]); 
} 
 
//----------------------------------------------- 
MultiRange::size_type MultiRange::size () const 
{ 
  return (m_ranges.size ()); 
} 
 
MultiRange::size_type MultiRange::cardinality () const 
{ 
  size_type result = 0; 
 
  for (size_type i = 0; i < m_ranges.size (); ++i) 
    { 
      const Range& r = m_ranges[i]; 
 
      result += r.cardinality (); 
    } 
 
  return (result); 
} 
 
MultiRange::size_type MultiRange::cardinalityUpTo (const ExpandedIdentifier& id) const
{
    // Loop over ranges in MultiRange and calculate hash for each
    // range

    size_type result = 0;
    for (unsigned int i = 0; i < m_ranges.size(); ++i) {
	const Range& range = m_ranges[i];
	result += range.cardinalityUpTo(id);
    }
    return (result);
}


//----------------------------------------------- 
bool MultiRange::has_overlap () const 
{ 
  range_vector::size_type i; 
  range_vector::size_type j; 
 
  for (i = 0; i < m_ranges.size (); ++i) 
    { 
      const Range& r = m_ranges[i]; 
      for (j = i + 1; j < m_ranges.size (); ++j) 
        { 
          const Range& s = m_ranges[j]; 
          if (r.overlaps_with (s)) return (true); 
        } 
    } 
 
  return (false); 
} 
 
//----------------------------------------------- 
MultiRange::identifier_factory MultiRange::factory_begin (bool sort) 
{ 
  const MultiRange& me = *this; 
  if (sort) m_it_count++;
  return (identifier_factory (me, sort)); 
} 
 
//----------------------------------------------- 
MultiRange::const_identifier_factory MultiRange::factory_begin (bool sort) const 
{ 
  const MultiRange& me = *this; 
  if (sort) m_it_count++;
  return (const_identifier_factory (me, sort)); 
} 
 
//----------------------------------------------- 
MultiRange::identifier_factory MultiRange::factory_end () 
{ 
  static identifier_factory factory;
 
  return (factory); 
} 
 
//----------------------------------------------- 
MultiRange::const_identifier_factory MultiRange::factory_end () const 
{ 
  static const_identifier_factory factory; 
 
  return (factory); 
} 
 
//----------------------------------------------- 
MultiRange::identifier_factory::identifier_factory () 
    : 
    m_sort(false),
    m_multirange(0)
{ 
}

//----------------------------------------------- 
MultiRange::identifier_factory::identifier_factory (const identifier_factory& other) 
    :
    m_id 		(other.m_id),
    m_sort 		(other.m_sort),
    m_id_fac_it  	(other.m_id_fac_it),
    m_id_fac_end 	(other.m_id_fac_end),
    m_range_it 		(other.m_range_it),
    m_range_end		(other.m_range_end),
    m_id_vec_it 	(other.m_id_vec_it),
    m_id_vec_end 	(other.m_id_vec_end),
    m_multirange	(other.m_multirange)
{
	    
    if (m_sort && m_multirange) {
	m_multirange->m_it_count++;
//  	std::cout << "MultiRange::identifier_factory::operator =: it_count" 
//  		  << m_multirange->m_it_count 
//  		  << std::endl;
    }
} 
 
//----------------------------------------------- 
MultiRange::identifier_factory::identifier_factory (const MultiRange& multirange, bool sort) 
    :
    m_sort(sort),
    m_range_it(multirange.m_ranges.begin()),
    m_range_end(multirange.m_ranges.end()),
    m_multirange(&multirange)
{ 
    if (m_range_it == m_range_end)return;  // no ranges
    /** 
     *  For sorted ids, must set up initial id vector. For unsorted,
     *  one can just set up iterators over ranges and ids.
     */
    if (m_sort) {
	if (m_multirange->m_ids.size() > 0) {
	    // id vector has already been created
	    // Set up iterator over ids
	    m_id_vec_it  = m_multirange->m_ids.begin();
	    m_id_vec_end = m_multirange->m_ids.end();
	    if (m_id_vec_it != m_id_vec_end) {
		// Set id
		m_id = *m_id_vec_it;
	    }
	}
	else {
	    // Limit sorted access to just under 500K ids
	    size_type nids = m_multirange->cardinality();
	    if (500000 < nids) {
		std::cout << "MultiRange::identifier_factory: unable to provide " 
			  << " access to sorted identifier iterator. Limited to < 500k ids." 
			  << std::endl;
	    }
	    else {
		// Access and sort ids, removing duplicates
		std::set<ExpandedIdentifier> ids;
		for (; m_range_it != m_range_end; ++m_range_it) {
		    m_id_fac_it  = (*m_range_it).factory_begin();
		    m_id_fac_end = (*m_range_it).factory_end();
		    for (; m_id_fac_it != m_id_fac_end; ++m_id_fac_it) {
			ids.insert(*m_id_fac_it);
		    }
		}
		std::set<ExpandedIdentifier>::iterator    firstIds = ids.begin();
		std::set<ExpandedIdentifier>::iterator    lastIds  = ids.end();
		for(; firstIds != lastIds; ++firstIds) {
		    m_multirange->m_ids.push_back (*firstIds); 
		}
		// Set up iterator over ids
		m_id_vec_it  = m_multirange->m_ids.begin();
		m_id_vec_end = m_multirange->m_ids.end();
		if (m_id_vec_it != m_id_vec_end) {
		    // Set id
		    m_id = *m_id_vec_it;
		}
	    }	    
	}
    }
    else {
	/** 
	 *  No sort: Set up the iterators for range and identifiers
	 */
	if (m_range_it != m_range_end) {
	    m_id_fac_it  = (*m_range_it).factory_begin();
	    m_id_fac_end = (*m_range_it).factory_end();
	    if(m_id_fac_it != m_id_fac_end) {
		// Set id
		m_id = *m_id_fac_it;
	    }
	}
    }
}
 
MultiRange::identifier_factory::~identifier_factory ()
{
    // remove id buffer in multirange if no iterator is still using it
    if (m_sort && m_multirange) {
	if(m_multirange->m_it_count) {
	    m_multirange->m_it_count--;
	}
	else {
	    std::cout << "MultiRange::identifier_factory::~identifier_factory: it_count ALREADY 0" 
		      << std::endl;
	}
//  	std::cout << "MultiRange::identifier_factory::~identifier_factory: it_count" 
//  		  << m_multirange->m_it_count 
//  		  << std::endl;
	if (!m_multirange->m_it_count) {
	    m_multirange->m_ids.clear();
//  	std::cout << "MultiRange::identifier_factory::~identifier_factory: cleared id buffer" 
//  		  << std::endl;
	}    
    }
}


//----------------------------------------------- 
MultiRange::identifier_factory& MultiRange::identifier_factory::operator = (const identifier_factory& other) 
{ 
  if (this != &other) {
    m_id 		= other.m_id;
    m_sort 		= other.m_sort;
    m_id_fac_it  	= other.m_id_fac_it;
    m_id_fac_end 	= other.m_id_fac_end;
    m_range_it 		= other.m_range_it;
    m_range_end		= other.m_range_end;
    m_id_vec_it 	= other.m_id_vec_it;
    m_id_vec_end 	= other.m_id_vec_end;
    m_multirange	= other.m_multirange;
    // Must increment number of iterators looking at the id buffer in
    // the multirange
    if (m_sort && m_multirange) {
	m_multirange->m_it_count++;
//  	std::cout << "MultiRange::identifier_factory::operator =: it_count" 
//  		  << m_multirange->m_it_count 
//  		  << std::endl;
    }
  }
  return (*this); 
} 
 
//----------------------------------------------- 
void MultiRange::identifier_factory::operator ++ () 
{ 
 
    if (m_id.fields () == 0) return; 

    /** 
     *  Advance iterators depending upon whether or not we have sorted.
     */

    if (m_sort) {
	m_id.clear();
	if (m_id_vec_it != m_id_vec_end) {
	    ++m_id_vec_it;
	    if (m_id_vec_it != m_id_vec_end) {
		m_id = *m_id_vec_it;
	    }
	}
    }
    else {
	m_id.clear();
	if (m_range_it != m_range_end) {
	    if (m_id_fac_it != m_id_fac_end) {
		++m_id_fac_it;
	    }
	    if (m_id_fac_it == m_id_fac_end) {
		++m_range_it;
		if (m_range_it != m_range_end) {
		    m_id_fac_it  = (*m_range_it).factory_begin();
		    m_id_fac_end = (*m_range_it).factory_end();
		}
	    }
	    if (m_id_fac_it != m_id_fac_end) {
		m_id = *m_id_fac_it;
	    }
	}
    }
} 

 
//----------------------------------------------- 
const ExpandedIdentifier& MultiRange::identifier_factory::operator * () const 
{ 
  return (m_id); 
} 
 
//----------------------------------------------- 
bool MultiRange::identifier_factory::operator == (const identifier_factory& other) const 
{ 
  if (m_id == other.m_id) return (true); 
  return (false); 
} 
 
//----------------------------------------------- 
bool MultiRange::identifier_factory::operator != (const identifier_factory& other) const 
{ 
  if (m_id != other.m_id) return (true); 
  return (false); 
} 

//----------------------------------------------- 
MultiRange::const_identifier_factory::const_identifier_factory () 
    : 
    m_sort(false),
    m_multirange(0)
{ 
} 

//----------------------------------------------- 
MultiRange::const_identifier_factory::const_identifier_factory (const const_identifier_factory& other) 
    :
    m_id 		(other.m_id),
    m_sort 		(other.m_sort),
    m_id_fac_it  	(other.m_id_fac_it),
    m_id_fac_end 	(other.m_id_fac_end),
    m_range_it 		(other.m_range_it),
    m_range_end		(other.m_range_end),
    m_id_vec_it 	(other.m_id_vec_it),
    m_id_vec_end 	(other.m_id_vec_end),
    m_multirange	(other.m_multirange)
{
	    
    if (m_sort && m_multirange) {
	m_multirange->m_it_count++;
//  	std::cout << "MultiRange::const_identifier_factory::operator =: it_count" 
//  		  << m_multirange->m_it_count 
//  		  << std::endl;
    }
} 
 
 
//----------------------------------------------- 
MultiRange::const_identifier_factory::const_identifier_factory (const MultiRange& multirange, bool sort) 
    :
    m_sort(sort),
    m_range_it(multirange.m_ranges.begin()),
    m_range_end(multirange.m_ranges.end()),
    m_multirange(&multirange)
{ 

//      std::cout << "MultiRange::const_identifier_factory::const_identifier_factory: sort "  
//  	      << sort << " size " << m_multirange->m_ids.size()
//  	      << " range it " << (m_range_it == m_range_end)
//  	      << " multi range " << (std::string)(*m_multirange)
//  	      << std::endl;
    

    if (m_range_it == m_range_end)return;  // no ranges
    /** 
     *  For sorted ids, must set up initial id vector. For unsorted,
     *  one can just set up iterators over ranges and ids.
     */
    if (m_sort) {
	if (m_multirange->m_ids.size() > 0) {
	    // id vector has already been created
	    // Set up iterator over ids
	    m_id_vec_it  = m_multirange->m_ids.begin();
	    m_id_vec_end = m_multirange->m_ids.end();
	    if (m_id_vec_it != m_id_vec_end) {
		// Set id
		m_id = *m_id_vec_it;
	    }
	}
	else {
	    // Limit sorted access to just under 500K ids
	    size_type nids = m_multirange->cardinality();
	    if (500000 < nids) {
		std::cout << "MultiRange::identifier_factory: unable to provide " 
			  << " access to sorted identifier iterator. Limited to < 500k ids." 
			  << std::endl;
	    }
	    else {
		// Access and sort ids, removing duplicates
		std::set<ExpandedIdentifier> ids;
		for (; m_range_it != m_range_end; ++m_range_it) {
		    m_id_fac_it  = (*m_range_it).factory_begin();
		    m_id_fac_end = (*m_range_it).factory_end();
		    for (; m_id_fac_it != m_id_fac_end; ++m_id_fac_it) {
			ids.insert(*m_id_fac_it);
		    }
		}
		std::set<ExpandedIdentifier>::iterator    firstIds = ids.begin();
		std::set<ExpandedIdentifier>::iterator    lastIds  = ids.end();
		for(; firstIds != lastIds; ++firstIds) {
		    m_multirange->m_ids.push_back (*firstIds); 
		}
		// Set up iterator over ids
		m_id_vec_it  = m_multirange->m_ids.begin();
		m_id_vec_end = m_multirange->m_ids.end();
		if (m_id_vec_it != m_id_vec_end) {
		    // Set id
		    m_id = *m_id_vec_it;
		}
	    }	    
	}
    }
    else {
	/** 
	 *  No sort: Set up the iterators for range and identifiers
	 */
	if (m_range_it != m_range_end) {
	    m_id_fac_it  = (*m_range_it).factory_begin();
	    m_id_fac_end = (*m_range_it).factory_end();
	    if(m_id_fac_it != m_id_fac_end) {
		// Set id
		m_id = *m_id_fac_it;
	    }
	}
    }
} 

 
MultiRange::const_identifier_factory::~const_identifier_factory ()
{
    // remove id buffer in multirange if no iterator is still using it
    if (m_sort && m_multirange) {
	if(m_multirange->m_it_count) {
	    m_multirange->m_it_count--;
	}
	else {
	    std::cout << "MultiRange::const_identifier_factory::~const_identifier_factory: it_count ALREADY 0" 
		      << std::endl;
	}
//  	std::cout << "MultiRange::const_identifier_factory::~const_identifier_factory: it_count" 
//  		  << m_multirange->m_it_count 
//  		  << std::endl;
	if (!m_multirange->m_it_count) {
	    m_multirange->m_ids.clear();
//  	std::cout << "MultiRange::const_identifier_factory::const_~identifier_factory: cleared id buffer" 
//  		  << std::endl;
	}    
    }
}

 
//----------------------------------------------- 
MultiRange::const_identifier_factory& MultiRange::const_identifier_factory::operator = (const const_identifier_factory& other) 
{ 
  if (this != &other) {
    m_id 		= other.m_id;
    m_sort 		= other.m_sort;
    m_id_fac_it  	= other.m_id_fac_it;
    m_id_fac_end 	= other.m_id_fac_end;
    m_range_it 		= other.m_range_it;
    m_range_end		= other.m_range_end;
    m_id_vec_it 	= other.m_id_vec_it;
    m_id_vec_end 	= other.m_id_vec_end;
    m_multirange	= other.m_multirange;
    // Must increment number of iterators looking at the id buffer in
    // the multirange
    if (m_sort && m_multirange) {
	m_multirange->m_it_count++;
//  	std::cout << "MultiRange::const_identifier_factory::operator =: it_count" 
//  		  << m_multirange->m_it_count 
//  		  << std::endl;
    }
  } 
  return (*this); 
} 
 
//----------------------------------------------- 
void MultiRange::const_identifier_factory::operator ++ () 
{ 
 
    if (m_id.fields () == 0) return; 

    /** 
     *  Advance iterators depending upon whether or not we have sorted.
     */

    if (m_sort) {
	m_id.clear();
	if (m_id_vec_it != m_id_vec_end) {
	    ++m_id_vec_it;
	    if (m_id_vec_it != m_id_vec_end) {
		m_id = *m_id_vec_it;
	    }
	}
    }
    else {
	m_id.clear();
	if (m_range_it != m_range_end) {
	    if (m_id_fac_it != m_id_fac_end) {
		++m_id_fac_it;
	    }
	    if (m_id_fac_it == m_id_fac_end) {
		++m_range_it;
		if (m_range_it != m_range_end) {
		    m_id_fac_it  = (*m_range_it).factory_begin();
		    m_id_fac_end = (*m_range_it).factory_end();
		}
	    }
	    if (m_id_fac_it != m_id_fac_end) {
		m_id = *m_id_fac_it;
	    }
	}
    }
} 
 
//----------------------------------------------- 
const ExpandedIdentifier& MultiRange::const_identifier_factory::operator * () const 
{ 
  return (m_id); 
} 
 
//----------------------------------------------- 
bool MultiRange::const_identifier_factory::operator == (const const_identifier_factory& other) const 
{ 
  if (m_id == other.m_id) return (true); 
  return (false); 
} 
 
//----------------------------------------------- 
bool MultiRange::const_identifier_factory::operator != (const const_identifier_factory& other) const 
{ 
  if (m_id != other.m_id) return (true); 
  return (false); 
} 
 

/** 
 *   This is a debugging facility: 
 *   two vectors of ExpandedIdentifiers are filled in : 
 *     one with sorted and cleaned from duplicates 
 *     the other contains only duplicates (empty when there is no overlap) 
 * 
 *   This function may be used for validity checks 
 */ 
void MultiRange::show_all_ids (std::vector <ExpandedIdentifier>& unique_ids, 
                               std::vector <ExpandedIdentifier>& duplicate_ids) const 
{ 
  range_vector::const_iterator i; 
 
  std::vector <ExpandedIdentifier> ids; 
 
  for (i = m_ranges.begin (); i != m_ranges.end (); ++i) 
    { 
      const Range& r = *i; 
      Range::const_identifier_factory f = r.factory_begin (); 
      while (f != r.factory_end ()) 
        { 
          const ExpandedIdentifier id = *f; 
 
          ids.push_back (id); 
 
          ++f; 
        } 
    } 
 
  std::sort (ids.begin (), ids.end ()); 
 
  std::vector<ExpandedIdentifier>::const_iterator it; 
 
  bool first = true; 
  ExpandedIdentifier previous; 
 
  for (it = ids.begin (); it < ids.end (); ++it) 
    { 
      const ExpandedIdentifier id = *it; 
 
      if (first) 
        { 
          first = false; 
          unique_ids.push_back (id); 
          previous = id; 
        } 
      else 
        { 
          if (id.match (previous)) 
            { 
              duplicate_ids.push_back (id); 
            } 
          else 
            { 
              unique_ids.push_back (id); 
              previous = id; 
            } 
        } 
    } 
 
  std::cout << unique_ids.size () << " unique ids " << 
    duplicate_ids.size () << " duplicate ids" << std::endl; 
} 
 
/** 
 *    The reduction algorithm is no longer provided, since it does not 
 *   work in the most general case. 
 *    See previous and temptative implementation in the RangeReduction.cxx  
 *   source file. 
 */ 
void MultiRange::reduce () 
{ 
} 
 
//----------------------------------------------- 
void MultiRange::show () const 
{
  show (std::cout);
}

void MultiRange::show (std::ostream& s) const 
{ 
  range_vector::size_type i; 
 
  for (i = 0; i < m_ranges.size (); ++i) 
    { 
      if (i > 0) s << std::endl; 
 
      const Range& r = m_ranges[i]; 
      r.show (s); 
    } 
} 
 
//----------------------------------------------- 
MultiRange::operator std::string () const 
{ 
  std::string result; 
 
  range_vector::size_type i; 
 
  for (i = 0; i < m_ranges.size (); ++i) 
    { 
      if (i > 0) result += " | "; 
 
      const Range& r = m_ranges[i]; 
      result += (std::string) r; 
    } 
 
  return (result); 
} 

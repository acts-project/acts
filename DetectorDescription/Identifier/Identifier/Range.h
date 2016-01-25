#ifndef __Range_h__ 
#define __Range_h__ 
 
#include <Identifier/ExpandedIdentifier.h> 
#include <cassert>
 
/** 
 *    A Range describes the possible ranges for the field values of an ExpandedIdentifier 
 * 
 *    Specifications can be : 
 *       No bound      * 
 *       Low bound     n: 
 *       High bound    :m 
 *       Both bounds   n:m 
 *       Enumeration   v1, v2, v3, ... , vn 
 * 
 *    Trailing * are implicit for all trailing fields 
 * 
 */ 
class Range 
{ 
public: 
 
  typedef ExpandedIdentifier::element_type element_type; 
  typedef ExpandedIdentifier::size_type size_type; 
 
  /** 
   *   This is the individual specification for the range of one ExpandedIdentifier field.  
   */ 
  class field 
  { 
    public : 
 
    /** 
     *   Characterizes the four possible modes of any field specification 
     */ 
    typedef enum 
    { 
      unbounded, 
      low_bounded, 
      high_bounded, 
      both_bounded, 
      enumerated 
    } mode; 
 
    typedef enum 
    { 
      none, 
      has_next, 
      has_previous, 
      has_both,
      has_wrap_around
    } continuation_mode; 
 
    typedef std::vector <element_type> element_vector; 
    typedef std::vector <size_type>    index_vector; 
 
    /// Create a wild-card value. 
    field (); 
 
    /// Create a field copy 
    field (const field& other); 

    /// Move constructor.
    field (field&& other);
 
    /// Create a unique value (understood as : low bound = high bound = value) 
    field (element_type value); 
 
    /// Create a full range specification (with explicit min and max) 
    field (element_type minimum, element_type maximum); 
 
    /// Some combined query functions on the specification mode 
    bool is_valued () const; 
    bool has_minimum () const; 
    bool has_maximum () const; 
    bool wrap_around () const; 
 
    /// Query the values 
    mode get_mode () const; 
    element_type get_minimum () const; 
    element_type get_maximum () const; 
    const element_vector& get_values () const; 
    ///  Returns false if previous/next is at end of range, or not possible
    bool get_previous (element_type current, element_type& previous) const; 
    bool get_next     (element_type current, element_type& next) const; 
    size_type get_indices () const; 
    index_vector get_indexes () const;      
    size_type get_bits () const; 
    element_type get_value_at (size_type index) const; 
    size_type get_value_index (element_type value) const; 
 
    /// Check if this is a pure wild card field 
    bool match_any () const; 
 
    /// The basic match operation 
    bool match (element_type value) const; 
 
    /// Check whether two fields overlap 
    bool overlaps_with (const field& other) const; 
 
    /// Set methods 
    void clear (); 
    void set (element_type minimum, element_type maximum); 
    void set_minimum (element_type value); 
    void set_maximum (element_type value); 
    void add_value (element_type value); 
    void set (const element_vector& values); 
    void set (bool wraparound); 
    void set_next (int next);
    void set_previous (int previous);
    field& operator = (const field& other); 
    void operator |= (const field& other); 
 
    operator std::string () const; 
    bool operator == (const field& other) const; 
    bool operator != (const field& other) const; 

    void show() const;
      
    /// Optimize - try to switch mode to both_bounded, set up lookup
    /// table for finding index from value
    void optimize();

  private : 
 
    typedef enum 
    { 
      max_indexes = 100
    } max_values; 

    /// Check mode - switch from enumerated to both_bounded if possible
    void check_for_both_bounded();

    /// Create index table from value table
    void create_index_table();

    /// Set m_indices
    void set_indices();

    element_type m_minimum; 
    element_type m_maximum; 
    element_vector m_values; 
    index_vector m_indexes; 
    size_type    m_indices;
    element_type m_previous; 
    element_type m_next; 
    mode m_mode; 
    continuation_mode m_continuation_mode; 
  }; 
 
  typedef std::vector<field> field_vector; 
 
  /** 
   *    This factory is able to generate all possible identifiers, from a  
   *  fully bounded Range. 
   *    The precondition is that the Range used to parameterize the factory  
   *  must have all its fields completely bounded. 
   */ 
  class identifier_factory 
  { 
  public: 
    identifier_factory (); 
    identifier_factory (const Range& range); 
 
    identifier_factory& operator = (const identifier_factory& other); 
 
    void operator ++ (); 
 
    const ExpandedIdentifier& operator * () const; 
    bool operator == (const identifier_factory& other) const; 
    bool operator != (const identifier_factory& other) const; 
 
  private: 
    std::vector<size_type> m_indices; 
    ExpandedIdentifier m_id; 
    ExpandedIdentifier m_min; 
    ExpandedIdentifier m_max; 
    const Range* m_range; 
  }; 
 
  class const_identifier_factory 
  { 
  public: 
    const_identifier_factory (); 
    const_identifier_factory (const Range& range); 
 
    const_identifier_factory& operator = (const const_identifier_factory& other); 
 
    void operator ++ (); 
 
    const ExpandedIdentifier& operator * () const; 
    bool operator == (const const_identifier_factory& other) const; 
    bool operator != (const const_identifier_factory& other) const; 
 
  private: 
    std::vector<size_type> m_indices; 
    ExpandedIdentifier m_id; 
    ExpandedIdentifier m_min; 
    ExpandedIdentifier m_max; 
    const Range* m_range; 
  }; 
 
  /// Constructors 
  Range (); 
  Range (const Range& other); 

  /// Assignment.
  Range& operator= (const Range& other);

  /// Constructor with move semantics.
  /// FIXME: Replace with rvalue reference once we can use C++11.
  Range (Range& other, bool); 
 
  /** 
   *   This is a sub-range copy constructor.  
   * It copies the portion of the other Range, starting from the  
   * specified starting index up to its last field. 
   */ 
  Range (const Range& other, size_type start); 
 
  /** 
   * Constructor with setup from a textual description. 
   * 
   *  The syntax is : 
   * 
   * range : 
   *      <value-range> [ "/" <value-range> ... ] 
   * 
   * value-range : 
   *      "*" 
   *    | <value> 
   *    | ":" <max> 
   *    | <min> ":" 
   *    | <min> ":" <max> 
   *    | <value> "," <value> "," ... "," <value> 
   * 
   */ 
//  explicit Range (const std::string& text); 
 
  /** 
   *   Construct from a simple ExpandedIdentifier. This implies that all fields 
   *   will have their min=max=id[i] 
   */ 
  Range (const ExpandedIdentifier& root); 
     
  /** 
   * Build a range from a textual description. 
   */ 
  void build (const std::string& text); 
 
  /** 
   *   Build a range from a single ExpandedIdentifier 
   *   (see similar constructor for comment) 
   */ 
  void build (const ExpandedIdentifier& root); 
 
  /// Modifications 
 
  void clear (); 
 
  /// Add a wild card field. 
  void add (); 
 
  /// Add a required value. (ie. low = high = value) 
  void add (element_type value); 
 
  /// Add a bounded value. 
  void add (element_type minimum, element_type maximum); 
 
  /// Add a range bounded by a minimum. 
  void add_minimum (element_type minimum); 
 
  /// Add a range bounded by a maximum. 
  void add_maximum (element_type maximum); 
 
  /// Add a range specified using a field  
  void add (const field& f); 
 
  /// Add a range specified using a field, with move semantics.
  void add (field&& f);
 
  /// Append a subrange 
  void add (const Range& subrange); 

  /// Append a subrange, with move semantics.
  void add (Range&& subrange); 
 
  /// Match an identifier 
  int match (const ExpandedIdentifier& id) const; 
 
  /// Accessors 
 
  /// Access the field elements 
  const field& operator [] (size_type index) const; 
  size_type fields () const; 
  bool is_empty () const; 
 
  /** 
   *   min and max ExpandedIdentifiers  
   *  (if they exist, ie. for fully bounded Ranges) 
   *  Question : what if the Range has wild cards ?? 
   */ 
  ExpandedIdentifier minimum () const; 
  ExpandedIdentifier maximum () const; 
 
  /** 
   *  Computes a possible cardinality : 
   *   - all bounded fields are counted as they are 
   *   - unbounded fields are conted for one value. 
   */ 
  size_type cardinality () const;
  //  Up to a given id
  size_type cardinalityUpTo (const ExpandedIdentifier& id) const;
  size_type cardinalityUpTo (const int* id) const;
 
  /// Identifier_factory management 
  identifier_factory factory_begin (); 
  const_identifier_factory factory_begin () const; 
  identifier_factory factory_end (); 
  const_identifier_factory factory_end () const; 
 
  /// Check if two Ranges overlap. 
  bool overlaps_with (const Range& other) const; 
 
  void show () const; 
  void show (std::ostream& s) const; 
 
  /// Produce a textual representation of the range using the input format 
  operator std::string () const; 
 
  bool operator == (const Range& other) const; 
  bool operator != (const Range& other) const; 

private: 
  field_vector m_fields; 
}; 
 
/** 
 *   A MultiRange combines several Ranges 
 */ 
class MultiRange 
{ 
public: 
 
  typedef std::vector<Range> range_vector; 
  typedef ExpandedIdentifier::element_type element_type; 
  typedef ExpandedIdentifier::size_type size_type; 
 
    /** 
   *    This factory is able to generate all possible identifiers, from a  
   *  fully bounded Range. 
   *    The precondition is that the Range used to parameterize the factory  
   *  must have all its fields completely bounded. 
   */ 
    class identifier_factory 
    { 
    public: 
	identifier_factory (); 
	identifier_factory (const identifier_factory& other); 
	identifier_factory (const MultiRange& multirange, bool sort); 
	~identifier_factory (); 
 
	identifier_factory& operator = (const identifier_factory& other); 

	void operator ++ (); 
 
	const ExpandedIdentifier& operator * () const; 
	bool operator == (const identifier_factory& other) const; 
	bool operator != (const identifier_factory& other) const; 
 
    private: 
 
	typedef std::vector<ExpandedIdentifier>	id_vec;
	typedef id_vec::iterator 		id_iterator;
	typedef id_vec::const_iterator 		id_const_iterator;

	ExpandedIdentifier		m_id;
	bool 				m_sort;
	Range::const_identifier_factory	m_id_fac_it; 
	Range::const_identifier_factory	m_id_fac_end; 
	range_vector::const_iterator   	m_range_it; 
	range_vector::const_iterator   	m_range_end; 
	id_iterator			m_id_vec_it;
	id_iterator			m_id_vec_end;
	const MultiRange*		m_multirange;
    }; 
 
    class const_identifier_factory 
    {
    public: 
	const_identifier_factory (); 
	const_identifier_factory (const const_identifier_factory& other); 
	const_identifier_factory (const MultiRange& multirange, bool sort); 
	~const_identifier_factory (); 
 
	const_identifier_factory& operator = (const const_identifier_factory& other); 
 
	void operator ++ (); 
 
	const ExpandedIdentifier& operator * () const; 
	bool operator == (const const_identifier_factory& other) const; 
	bool operator != (const const_identifier_factory& other) const; 
 
    private: 
	typedef std::vector<ExpandedIdentifier> id_vec;
	typedef id_vec::iterator 		id_iterator;
	typedef id_vec::const_iterator 		id_const_iterator;

	ExpandedIdentifier		m_id;
	bool 				m_sort;
	Range::const_identifier_factory	m_id_fac_it; 
	Range::const_identifier_factory	m_id_fac_end; 
	range_vector::const_iterator   	m_range_it; 
	range_vector::const_iterator   	m_range_end; 
	id_iterator			m_id_vec_it;
	id_iterator			m_id_vec_end;
	const MultiRange*		m_multirange;
    }; 
 
  /// Constructors 
  MultiRange (); 
  MultiRange (const MultiRange& other); 
  MultiRange (const std::string& text); 


  /// Assignment.
  MultiRange& operator= (const MultiRange& other);

  /** 
   *   Construct a non-overlapping MultiRange from 
   *   two overlapping ones 
   */ 
  MultiRange (const Range& r, const Range& s); 
 
  /// Build a range from a textual description. 
  void build (const std::string& text); 
 
  /// Modifications 
 
  void clear (); 
 
  void add (const Range& range); 

  /// Add with move semantics.
  void add (Range& range, bool); 

  /// Add a Range made from a single ExpandedIdentifier 
  void add (const ExpandedIdentifier& id); 

  /// Remove a Range made from a single ExpandedIdentifier
  void remove_range (const ExpandedIdentifier& id); 
 
  /// Create a new empty Range that can be adapted afterwards 
  Range& add_range (); 
 
  /// Get the last entered Range 
  Range& back (); 
 
  /// Match an identifier 
  int match (const ExpandedIdentifier& id) const; 
 
  /// Accessors 
  const Range& operator [] (size_type index) const; 
  size_type size () const; 
 
  /** 
   *  Computes a possible cardinality from all ranges. 
   */ 
  size_type cardinality () const; 
  //  Up to a given id
  size_type cardinalityUpTo (const ExpandedIdentifier& id) const;
 
  /// Check if there are overlaps between any couple of Ranges 
  bool has_overlap () const; 
  void reduce (); 
 
  // identifier_factory management 
  identifier_factory 		factory_begin (bool sort = false); 
  const_identifier_factory 	factory_begin (bool sort = false) const; 
  identifier_factory 		factory_end (); 
  const_identifier_factory 	factory_end () const; 

 
  void show () const; 
  void show (std::ostream& s) const; 
 
  /// Generate a textual representation of the multirange using the input format 
  operator std::string () const; 
 
  void show_all_ids (std::vector <ExpandedIdentifier>& unique_ids, 
                     std::vector <ExpandedIdentifier>& duplicate_ids) const; 
 
private: 
  friend class identifier_factory;
  friend class const_identifier_factory;
  typedef std::vector<ExpandedIdentifier>	id_vec;
  range_vector 		m_ranges; 
  // number of iterators accessing m_ids
  mutable size_type 	m_it_count;  
  mutable id_vec	m_ids;

}; 
 


//-------------------
// inline definitions
//-------------------


//------------------------------------------------------------------
inline Range::field::mode Range::field::get_mode () const 
//------------------------------------------------------------------
{ 
  return (m_mode); 
} 

//------------------------------------------------------------------
inline Range::element_type Range::field::get_minimum () const  
//------------------------------------------------------------------
{ 
  return (m_minimum); 
} 
 
//------------------------------------------------------------------
inline Range::element_type Range::field::get_maximum () const 
//------------------------------------------------------------------
{ 
  return (m_maximum); 
} 
 
//------------------------------------------------------------------
inline const Range::field::element_vector& Range::field::get_values () const 
//------------------------------------------------------------------
{ 
  return (m_values); 
} 

//------------------------------------------------------------------
inline ExpandedIdentifier::size_type Range::field::get_indices () const 
//------------------------------------------------------------------
{ 
    return (m_indices);
}

//------------------------------------------------------------------
inline Range::field::index_vector Range::field::get_indexes () const
//------------------------------------------------------------------
{
    return (m_indexes);
}


//------------------------------------------------------------------
inline ExpandedIdentifier::size_type Range::field::get_bits () const 
//------------------------------------------------------------------
{
  ExpandedIdentifier::size_type result = 1;

  size_t indices = get_indices ();

  indices--;
  if (indices > 0)
    {
      result = 0;
      while (indices > 0)
        {
          indices /= 2;
          result++;
        }
    }

  return (result);
}

//----------------------------------------------- 
inline Range::element_type 
Range::field::get_value_at (size_type index) const 
//----------------------------------------------- 
{ 
    // Only both_bounded and enumerated are valid to calculate the
    // value.
    // both_bounded if the more frequent case and so comes first.

    if (both_bounded == m_mode) {
	return (m_minimum + index); 
//  	if (index >= (size_type) (m_maximum - m_minimum + 1)) return (0); 
//  	else return (m_minimum + index); 
    }
    else if (enumerated == m_mode) {
	return (m_values[index]); 
//        if (index >= m_values.size ()) return (0); 
//        else return (m_values[index]); 
    }
 
    return (0); 
} 
 
//----------------------------------------------- 
inline ExpandedIdentifier::size_type 
Range::field::get_value_index (element_type value) const
//----------------------------------------------- 
{

    // Only both_bounded and enumerated are valid to calculate the
    // index.
    // both_bounded if the more frequent case and so comes first.

    if (both_bounded == m_mode) {
//  	if ((value >= m_minimum) &&
//  	    (value <= m_maximum)) {
	    return (value - m_minimum); 
//  	}
    }
    else if (enumerated == m_mode) {
//	if ((int)m_indexes.size() > ((int)value - (int)m_minimum)) {
	if (m_indexes.size()) {
	    // Table has been created, do simple lookup
            assert (value >= m_minimum && value - m_minimum < (int)m_indexes.size());
	    return (m_indexes[value - m_minimum]);
	}
	else {
	    for (size_type i = 0; i < m_values.size (); ++i) { 
		if (m_values[i] == value) return (i); 
	    }
	}
    }
 
    return (0); 
}

//----------------------------------------------- 
inline bool Range::field::match (element_type value) const 
//----------------------------------------------- 
{ 
    size_t i; 
 
    if (both_bounded == m_mode) {
	return ((value >= m_minimum) && 
		(value <= m_maximum)); 
    }
    else if (enumerated == m_mode) {
	for (i = 0; i < m_values.size (); ++i) { 
	    if (value == m_values[i]) return (true); 
	} 
	return (false); 
    }
    else if (unbounded == m_mode) {
	return (true); 
    }
    else if (high_bounded == m_mode) {
	return (value <= m_maximum); 
    }
    else if (low_bounded == m_mode) {
	return (value >= m_minimum); 
    }
    return (false);
} 

//----------------------------------------------- 
inline Range::size_type Range::fields () const 
//----------------------------------------------- 
{ 
  return (m_fields.size ()); 
} 
 
//--------------------------------------------------------------------------
inline const Range::field& Range::operator [] (Range::size_type index) const 
//--------------------------------------------------------------------------
{ 
  if (index >= m_fields.size ()) 
    { 
      static const field f; 
 
      return (f); 
    } 
 
  return (m_fields[index]); 
} 
 
//----------------------------------------------- 
inline bool Range::is_empty () const 
//----------------------------------------------- 
{ 
  if (m_fields.size () == 0) return (true); 
  return (false); 
} 


/**
 *   Get the cardinality from the beginning up to the given Identifier
 *   expanded into a int array. 
 */
//--------------------------------------------------------------------------
inline Range::size_type Range::cardinalityUpTo (const int* id) const 
//--------------------------------------------------------------------------
{ 
  size_type result = 0; 

  const Range& me = *this; 
  size_type level = 0;
  for (; level < fields (); ++level) {

      const field& f = me[level]; 

      size_type card = f.get_value_index (id[level]);
      
      for (size_type k = level + 1; k < fields(); ++k) {

	  const field& f = me[k]; 

	  card *= f.get_indices();
      }
      result += card;
  }
  return result;
} 


#endif 

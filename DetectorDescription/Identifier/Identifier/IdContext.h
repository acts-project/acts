/***************************************************************************
 Identifier package
 -----------------------------------------
 Copyright (C) 2002 by ATLAS Collaboration
 ***************************************************************************/

//<doc><file>	$Id: IdContext.h,v 1.3 2003-08-05 17:03:19 schaffer Exp $
//<version>	$Name: not supported by cvs2svn $

#ifndef IDENTIFIER_IDCONTEXT_H
# define IDENTIFIER_IDCONTEXT_H

//<<<<<< INCLUDES                                                       >>>>>>

#include "Identifier/ExpandedIdentifier.h"

//<<<<<< PUBLIC DEFINES                                                 >>>>>>
//<<<<<< PUBLIC CONSTANTS                                               >>>>>>
//<<<<<< PUBLIC TYPES                                                   >>>>>>
//<<<<<< PUBLIC VARIABLES                                               >>>>>>
//<<<<<< PUBLIC FUNCTIONS                                               >>>>>>
//<<<<<< CLASS DECLARATIONS                                             >>>>>>

/**
 *  class IdContext
 *
 *  This class saves the "context" of an expanded identifier
 *  (ExpandedIdentifier) for compact or hash versions (Identifier32 or
 *  IdentifierHash). This context is composed of 
 *
 *    1) begin and end indices of fields that are stored in the
 *       compact/hash id
 *    2) a possible "prefix" identifier for cases where the begin
 *       index is not 0 or the top level of the expaneded identifier. 
 *
 *  The IdContext is needed when only some of the identifier levels
 *   are to encoded in the compact/hash ids. 
 */

class IdContext
{
public:

    //
    // Define public typedefs
    //
    typedef ExpandedIdentifier::size_type 	size_type;

    // default constructor
    IdContext();
    // constructor with full initialization
    IdContext(const ExpandedIdentifier& prefix, 
	      size_type begin_index, 
	      size_type end_index);

    //
    // accessors
    //
    const ExpandedIdentifier&		prefix_id	(void) const;

    // indices of the first/last identifier fields
    size_type				begin_index	(void) const;
    size_type				end_index 	(void) const;
    
    //
    // modifiers
    //
    void				set	(const ExpandedIdentifier& prefix,
						 size_type begin_index,
						 size_type end_index);
    
private:
    
    ExpandedIdentifier	m_prefix;
    size_type		m_begin_index;
    size_type		m_end_index;
};

    


//<<<<<< INLINE PUBLIC FUNCTIONS                                        >>>>>>
//<<<<<< INLINE MEMBER FUNCTIONS                                        >>>>>>

inline IdContext::IdContext()
    :
    m_begin_index(0),
    m_end_index(0)
{}

inline IdContext::IdContext(const ExpandedIdentifier& prefix, 
			    size_type begin_index, 
			    size_type end_index)
    :
    m_prefix(prefix),
    m_begin_index(begin_index),
    m_end_index(end_index)
{}

inline const ExpandedIdentifier&		
IdContext::prefix_id	(void) const
{
    return (m_prefix);
}

inline IdContext::size_type			
IdContext::begin_index	(void) const
{
    return (m_begin_index);
}

inline IdContext::size_type			
IdContext::end_index 	(void) const
{
    return (m_end_index);
}

inline void			
IdContext::set (const ExpandedIdentifier& prefix,
		size_type begin_index,
		size_type end_index)
{
    m_prefix    	= prefix;
    m_begin_index 	= begin_index;
    m_end_index 	= end_index;
}

#endif // IDENTIFIER_IDCONTEXT_H

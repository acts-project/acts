/***************************************************************************
 Identifier package
 -----------------------------------------
 Copyright (C) 2002 by ATLAS Collaboration
 ***************************************************************************/

//<doc><file>	$Id: IdHelper.h,v 1.6 2004-04-21 09:57:48 fledroit Exp $
//<version>	$Name: not supported by cvs2svn $

#ifndef IDENTIFIER_IDCONVERSIONSTRATEGY_H
# define IDENTIFIER_IDCONVERSIONSTRATEGY_H

//<<<<<< INCLUDES                                                       >>>>>>
#include <string>                                   
//<<<<<< PUBLIC DEFINES                                                 >>>>>>
//<<<<<< PUBLIC CONSTANTS                                               >>>>>>
//<<<<<< PUBLIC TYPES                                                   >>>>>>

class Identifier;
class IdentifierHash;
class IdContext;
class IdDictMgr;
class IMessageSvc;

//<<<<<< PUBLIC VARIABLES                                               >>>>>>
//<<<<<< PUBLIC FUNCTIONS                                               >>>>>>
//<<<<<< CLASS DECLARATIONS                                             >>>>>>

//
//  class  IdHelper
//
//  This is an abstract base class for helper classes that know how to
//  convert Identifier <-> IdentifierHash objects.
//
class  IdHelper
{
public:
    
    virtual ~IdHelper(void);

    // Create compact id from hash id (return == 0 for OK)
    virtual int         get_id          (const IdentifierHash& hash_id,
                                         Identifier& id,
                                         const IdContext* context = 0) const = 0;
    
    // Create hash id from compact id (return == 0 for OK)
    virtual int         get_hash        (const Identifier& id, 
                                         IdentifierHash& hash_id,
                                         const IdContext* context = 0) const = 0;

    // Initialization from the identifier dictionary
    virtual int         initialize_from_dictionary(const IdDictMgr& dict_mgr) = 0;

    // retrieve version of the dictionary
    virtual std::string   dictionaryVersion  (void) const= 0;

    /// Checks are performed by default in debug compilation and NOT
    /// in optimized compilation. One can switch or query this mode for
    /// any idHelper with the following methods:
    virtual bool		do_checks	(void) const = 0;
    virtual void		set_do_checks	(bool do_checks) const = 0;
    /// Neighbour initialization is performed by default
    /// One can switch or query this mode for
    /// any idHelper with the following method:
    virtual bool		do_neighbours   	(void) const = 0;
    virtual void		set_do_neighbours	(bool do_neighbours) const = 0;

    // setting pointer to the MessageSvc
    virtual void                setMessageSvc  (IMessageSvc* msgSvc) = 0;

    virtual void                setDictVersion  (const IdDictMgr& dict_mgr, const std::string& name) = 0;
};





//<<<<<< INLINE PUBLIC FUNCTIONS                                        >>>>>>
//<<<<<< INLINE MEMBER FUNCTIONS                                        >>>>>>

#endif // IDENTIFIER_IDCONVERSIONSTRATEGY_H

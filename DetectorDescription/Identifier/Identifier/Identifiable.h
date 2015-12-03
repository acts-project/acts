/***************************************************************************
 Identifier Package
 -----------------------------------------
 Copyright (C) 1998 by ATLAS Collaboration
 ***************************************************************************/

//<doc><file>	$Id: Identifiable.h,v 1.4 2002-09-20 11:40:48 schaffer Exp $
//<version>	$Name: not supported by cvs2svn $

#ifndef IDENTIFIER_IDENTIFIABLE_H
# define IDENTIFIER_IDENTIFIABLE_H

//<<<<<< INCLUDES                                                       >>>>>>
//<<<<<< PUBLIC DEFINES                                                 >>>>>>
//<<<<<< PUBLIC CONSTANTS                                               >>>>>>
//<<<<<< PUBLIC TYPES                                                   >>>>>>

class Identifier;
class IdentifierHash;
class IdHelper;

//<<<<<< PUBLIC VARIABLES                                               >>>>>>
//<<<<<< PUBLIC FUNCTIONS                                               >>>>>>
//<<<<<< CLASS DECLARATIONS                                             >>>>>>

//
//  class Identifiable 
//
//  This class provides an abstract interface to an Identifiable
//  object.
//
//  It is "identifiable" in the sense that each object must have
//  an identify method returning an Identifier.
//  
//  The interface also is extended to also provide access to a "hash"
//  form of an identifier. And there is the possiblity to add a
//  conversion strategy to allow conversion from Identifier <->
//  IdentifierHash.
//

class Identifiable 
{
public:
    
    virtual ~Identifiable(void);
    
    virtual Identifier		identify() const = 0;

    virtual IdentifierHash	identifyHash() const;

    virtual const IdHelper* 	getHelper() const;
};


//<<<<<< INLINE PUBLIC FUNCTIONS                                        >>>>>>
//<<<<<< INLINE MEMBER FUNCTIONS                                        >>>>>>

#endif // IDENTIFIER_IDENTIFIABLE_H

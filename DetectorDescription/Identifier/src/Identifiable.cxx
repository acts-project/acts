/***************************************************************************
 Identifier Package
 -----------------------------------------
 Copyright (C) 2001 by ATLAS Collaboration
 ***************************************************************************/

//<doc><file>	$Id: Identifiable.cxx,v 1.4 2002-09-20 11:40:51 schaffer Exp $
//<version>	$Name: not supported by cvs2svn $

//<<<<<< INCLUDES                                                       >>>>>>

#include "Identifier/Identifiable.h"
#include "Identifier/IdentifierHash.h"
#include "Identifier/IdHelper.h"

Identifiable::~Identifiable ()
{
}

// default implementation
IdentifierHash	Identifiable::identifyHash() const
{
    IdentifierHash result;
    return (result);
}

// default implementation
const IdHelper* 
Identifiable::getHelper() const
{
    return (0);
}


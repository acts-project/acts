///////////////////////////////////////////////////////////////////
// MsgMacros.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_CORE_MSGMACROS_H
#define ACTS_CORE_MSGMACROS_H 1

#include "CoreInterfaces/MsgBase.h"
#define MSG_LVL_NOCHK(lvl, x)               \
        Acts::MsgBase::msg(lvl) << x << endmsg

#define MSG_LVL_CHK(lvl, x)                 \
        if (Acts::MsgBase::msgLvl(lvl)) Acts::MsgBase::msg(lvl) << x << endmsg

#define MSG_VERBOSE(x)  MSG_LVL_CHK(MSG::VERBOSE,x)
#define MSG_DEBUG(x)    MSG_LVL_CHK(MSG::DEBUG,x)
#define MSG_INFO(x)     MSG_LVL_NOCHK(MSG::INFO,x)  
#define MSG_WARNING(x)  MSG_LVL_NOCHK(MSG::WARNING,x)
#define MSG_ERROR(x)    MSG_LVL_NOCHK(MSG::ERROR,x) 
#define MSG_FATAL(x)    MSG_LVL_NOCHK(MSG::FATAL,x) 
#define MSG_ALWAYS(x)   MSG_LVL_NOCHK(MSG::ALWAYS,x)

#endif // ACTS_CORE_MSGMACROS_H

///////////////////////////////////////////////////////////////////
// MsgMacros.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_CORE_MSGMACROS_H
#define ATS_CORE_MSGMACROS_H 1

#ifdef ATS_CORE_MSG_PLUGIN
#include ATS_CORE_MSG_PLUGIN
#else 

#include "CoreInterfaces/MsgBase.h"
#define MSG_LVL_NOCHK(lvl, x)               \
        Ats::MsgBase::msg(lvl) << x << endmsg

#define MSG_LVL_CHK(lvl, x)                 \
        if (Ats::MsgBase::msgLvl(lvl)) Ats::MsgBase::msg(lvl) << x << endmsg

#define MSG_VERBOSE(x)  MSG_LVL_CHK(MSG::VERBOSE,x)
#define MSG_DEBUG(x)    MSG_LVL_CHK(MSG::DEBUG,x)
#define MSG_INFO(x)     MSG_LVL_NOCHK(MSG::INFO,x)  
#define MSG_WARNING(x)  MSG_LVL_NOCHK(MSG::WARNING,x)
#define MSG_ERROR(x)    MSG_LVL_NOCHK(MSG::ERROR,x) 
#define MSG_FATAL(x)    MSG_LVL_NOCHK(MSG::FATAL,x) 
#define MSG_ALWAYS(x)   MSG_LVL_NOCHK(MSG::ALWAYS,x)

#endif //ATS_GAUDI_BUILD

#endif // ATS_CORE_MSGMACROS_H

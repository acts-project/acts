///////////////////////////////////////////////////////////////////
// MsgMacros.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef BASECOMPSINTERFACES_MSGMACROS_H
#define BASECOMPSINTERFACES_MSGMACROS_H 1

// this is the ATHENA VERSION
#ifndef ATS_GAUDI_BUILD

// Athena version
#include "AthenaBaseComps/AthMsgStreamMacros.h"

#define MSG_VERBOSE(x)  ATH_MSG_VERBOSE(x)
#define MSG_DEBUG(x)    ATH_MSG_DEBUG(x) 
#define MSG_INFO(x)     ATH_MSG_INFO(x)   
#define MSG_WARNING(x)  ATH_MSG_WARNING(x)
#define MSG_ERROR(x)    ATH_MSG_ERROR(x)  
#define MSG_FATAL(x)    ATH_MSG_FATAL(x)  
#define MSG_ALWAYS(x)   ATH_MSG_ALWAYS(x) 

// GAUDI build
#else

#define MSG_LVL_NOCHK(lvl, x)               \
        this->msg(lvl) << x << endmsg

#define MSG_LVL_CHK(lvl, x)                 \
        if (this->msgLvl(lvl)) this->msg(lvl) << x << endmsg 

#define MSG_VERBOSE(x)  MSG_LVL_CHK(MSG::VERBOSE,x)
#define MSG_DEBUG(x)    MSG_LVL_CHK(MSG::DEBUG,x)
#define MSG_INFO(x)     MSG_LVL_NOCHK(MSG::INFO,x)  
#define MSG_WARNING(x)  MSG_LVL_NOCHK(MSG::WARNING,x)
#define MSG_ERROR(x)    MSG_LVL_NOCHK(MSG::ERROR,x) 
#define MSG_FATAL(x)    MSG_LVL_NOCHK(MSG::FATAL,x) 
#define MSG_ALWAYS(x)   MSG_LVL_NOCHK(MSG::ALWAYS,x)

#endif

// endif of file preprocessor include
#endif

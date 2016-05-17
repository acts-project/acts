///////////////////////////////////////////////////////////////////
// ExtrapolationCell.cpp, ACTS project 
///////////////////////////////////////////////////////////////////

#include "ACTS/Extrapolation/ExtrapolationCell.hpp"

std::vector<std::string> Acts::ExtrapolationCode::s_ecodeNames = { "Unset",                  
                                                                   "InProgress",
                                                                   "SuccessDestination",     
                                                                   "SuccessBoundaryReached", 
                                                                   "SuccessPathLimit", 
                                                                   "SuccessMaterialLimit",  
                                                                   "Recovered",             
                                                                   "FailureDestination",     
                                                                   "FailureLoop",     
                                                                   "FailureNavigation",      
                                                                   "FailureUpdateKill",      
                                                                   "FailureConfiguration",
                                                                   "LeftKnownWorld" };

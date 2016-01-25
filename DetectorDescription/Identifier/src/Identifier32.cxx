
#include "Identifier/Identifier32.h"
#include <stdarg.h>
#include <stdio.h>
#include <algorithm>

#include <iostream>
#include <iomanip>


//-----------------------------------------------
std::string Identifier32::getString() const
{
  std::string result;
  char temp[20];

  sprintf (temp, "0x%x", (unsigned int)m_id);
  result += temp;
  return (result);
}

//-----------------------------------------------
void Identifier32::show () const
{
    const Identifier32& me = *this;
    std::cout << me.getString();
}



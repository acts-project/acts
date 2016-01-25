
#include "Identifier/Identifier.h"
#include <stdarg.h>
#include <stdio.h>
#include <algorithm>

#include <iostream>
#include <iomanip>




//-----------------------------------------------
void Identifier::set (const std::string& id)
{
  sscanf (id.c_str(), "0x%" IDENTIFIER_PCODE "x", &m_id);
}


//-----------------------------------------------
std::string Identifier::getString() const
{
  std::string result;
  char temp[20];

  sprintf (temp, "0x%" IDENTIFIER_PCODE "x", (Identifier::value_type)m_id);
  result += temp;
  return (result);
}

//-----------------------------------------------
void Identifier::show () const
{
    const Identifier& me = *this;
    std::cout << me.getString();
}



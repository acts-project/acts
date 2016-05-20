// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_CHARGEDEFINITION_H
#define ACTS_CHARGEDEFINITION_H 1

namespace Acts
{
  class ChargedPolicy
  {
  public:
    ChargedPolicy(double charge):
      m_dCharge(charge){}

    ChargedPolicy(const ChargedPolicy& copy) = default;
    ChargedPolicy(ChargedPolicy&& moved) = default;
    ~ChargedPolicy() = default;
    ChargedPolicy& operator=(const ChargedPolicy& rhs) = default;
    ChargedPolicy& operator=(ChargedPolicy&& rhs) = default;

    bool operator==(const ChargedPolicy& rhs) const
    {
      return m_dCharge == rhs.m_dCharge;
    }

    bool operator!=(const ChargedPolicy& rhs) const
    {
      return !(*this == rhs);
    }

    double getCharge() const
    {
      return m_dCharge;
    }

    void setCharge(double charge)
    {
      m_dCharge = charge;
    }

    void flipSign()
    {
      m_dCharge *= -1.;
    }

  private:
    double m_dCharge;
  };

  class NeutralPolicy
  {
  public:
    bool operator==(const NeutralPolicy&) const
    {
      return true;
    }

    bool operator!=(const NeutralPolicy& rhs) const
    {
      return !(*this == rhs);
    }
    double getCharge() const
    {
      return 0.;
    }
  };
} // end of namespace Acts

#endif // ACTS_CHARGEDEFINITION_H

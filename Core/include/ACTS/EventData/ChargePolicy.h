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

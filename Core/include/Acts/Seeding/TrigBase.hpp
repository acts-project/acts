#include <cmath>


#define MAX_SILICON_LAYER_NUM   19
#define OffsetEndcapPixels  7
#define OffsetBarrelSCT     3
#define OffsetEndcapSCT    10

template <typename space_point_t> class TrigSiSpacePointBase {

 public:

  // using cylindrical co-ordinates, no. errors
  TrigSiSpacePointBase(long layer,
		       double r,  double phi,  double z,
		       double dr=0.0, double dz=0.0, const space_point_t* offlineSpacePoint = nullptr) :
    m_layer(layer),
    m_r(r), m_phi(phi), m_z(z),
    m_dr(dr), m_dz(dz), 
    m_offlineSpacePoint(offlineSpacePoint)
    {
      m_x = r * std::cos(phi);
      m_y = r * std::sin(phi);
      m_barCode=-1;
 


      if (m_offlineSpacePoint) {
        if (m_offlineSpacePoint->sourceLinks().size() == 1) {  // pixels have 1 SL
          m_isPixel = true;
        } 
        else {
          m_isPixel = false;
        }      
      }
      else {
        m_isPixel = false;//Arbitrarily choose value when no offline spacepoint
      }
    }

  // Destructor
    virtual ~TrigSiSpacePointBase() = default;
  
  // Methods to set data members
  void r(  const double r  ) {m_r   = r;  }
  void phi(const double phi) {m_phi = phi;}
  void z(  const double z  ) {m_z   = z;  }
  void x(  const double x  ) {m_x   = x;  }
  void y(  const double y  ) {m_y   = y;  }
  void dr(  const double dr  ) {m_dr   = dr;  }
  void dz(  const double dz  ) {m_dz   = dz;  }

  void barCode(int code) {m_barCode = code;}

  double r()    const {return m_r;}
  double phi()  const {return m_phi;}
  double z()    const {return m_z;}
  double dr()   const {return m_dr;} 
  double dz()   const {return m_dz;}
  double x()    const {return m_x;}
  double y()    const {return m_y;}
  long layer()  const {return m_layer;}

  bool isPixel() const {return  m_isPixel;}
  bool isSCT()   const {return !m_isPixel;}

  // Methods to calculate associated values

  double eta(double z0) const {
    double zr = (m_z-z0)/m_r; 
    return log(zr+std::sqrt(1.+zr*zr));
  }

  int barCode() const {return m_barCode;}
  const space_point_t* offlineSpacePoint() const {return m_offlineSpacePoint;}

 protected:

  long m_layer;

  double	m_r;
  double	m_x;
  double	m_y;
  
  double	m_phi;
  double	m_z;
  double	m_dr;
  double	m_dz;

  int           m_barCode;   //MC truth association
  bool          m_isPixel; //Stores whether spacepoint is Pixel or SCT

	const space_point_t* m_offlineSpacePoint;
};


template <typename space_point_t> class TrigInDetTriplet {

 public:
   TrigInDetTriplet() = delete; //to prevent creation w/o initialization

 TrigInDetTriplet(TrigSiSpacePointBase<space_point_t> s1, TrigSiSpacePointBase<space_point_t> s2, TrigSiSpacePointBase<space_point_t> s3, float Q) :
    m_s1(std::move(s1)), m_s2(std::move(s2)), m_s3(std::move(s3)), m_Q(Q) {};

 TrigInDetTriplet(TrigInDetTriplet* t) :
    m_s1(t->m_s1), m_s2(t->m_s2), m_s3(t->m_s3), m_Q(t->m_Q) {};

  const TrigSiSpacePointBase<space_point_t>& s1() const {return m_s1;}
  const TrigSiSpacePointBase<space_point_t>& s2() const {return m_s2;}
  const TrigSiSpacePointBase<space_point_t>& s3() const {return m_s3;}
  float Q() const {return m_Q;}
  void Q(double newQ) {m_Q = newQ;}

 protected:

  TrigSiSpacePointBase<space_point_t> m_s1;
  TrigSiSpacePointBase<space_point_t> m_s2;
  TrigSiSpacePointBase<space_point_t> m_s3;
  float m_Q;//Quality
};







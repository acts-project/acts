#include <stdint.h>
#include <map>
#include <ostream>
#include <iostream>
#include <atomic>
#include <vector> 

namespace Acts {

// class IRoiDescriptor {

//   public:

//     // iterator 
//     typedef std::vector<const IRoiDescriptor*>::const_iterator roi_iterator;

//   public:

//     /**
//     * @brief default constructor
//     */
//     IRoiDescriptor() { }
  
//     // Destructor
//     virtual ~IRoiDescriptor() { }

//     /// Methods to retrieve data members
    
//     /// directions
//     virtual   double phi()  const = 0; 
//     virtual   double eta()  const = 0;
//     virtual   double zed()  const = 0;

//     /// the zed and eta values at the most forward
//     /// and most rear ends of the RoI
//     virtual   double zedPlus()  const = 0;
//     virtual   double zedMinus() const = 0;

//     virtual   double etaPlus()  const = 0;
//     virtual   double etaMinus() const = 0;

//     /// extreme phi values
//     virtual   double phiPlus()  const = 0;
//     virtual   double phiMinus() const = 0;


//     /// identifiers
//     virtual   unsigned int roiId()   const = 0; 
//     virtual   unsigned int l1Id()    const = 0;
//     virtual   unsigned int roiWord() const = 0;


//     /// which roi version?
//     virtual int version() const = 0;

//     /// cast to a string
//     virtual operator std::string() const = 0;  


//     /// is this a full detector RoI?
//     virtual bool isFullscan() const = 0;


//     /// Super RoI access methods 

//     /// am I a SuperRoi?
//     virtual bool composite() const = 0;

//     /// number of constituents
//     virtual unsigned size() const = 0;

//     /// find an RoiDescriptor constituent
//     virtual const IRoiDescriptor* at(int i) const = 0;

//     /// const limit iterators
//     virtual roi_iterator  begin() const = 0;
//     virtual roi_iterator  end()   const = 0;


//     /// useful methods to determine whether items lie
//     /// partially within the RoI
    
//     /// accessors to calculate z position at radius r
//     /// along the RoI boundaries
//     virtual  double zedMin(double r) const = 0;
//     virtual  double zedMax(double r) const = 0;

//     /// accessors to calculate r position at position
//     /// z along the RoI boundaries
//     virtual  double rhoMin(double z) const = 0;
//     virtual  double rhoMax(double z) const = 0;


//     /// return the gradients 
//     virtual double dzdrMinus() const = 0;
//     virtual double dzdrPlus()  const = 0;

//     virtual double drdzMinus() const = 0;
//     virtual double drdzPlus()  const = 0;

//     /// zed limits at some outer radius
//     virtual   double zedOuterPlus()  const = 0;
//     virtual   double zedOuterMinus() const = 0;


// //   inline
// //   std::ostream& operator<<( std::ostream& s, const IRoiDescriptor& d ) {
// //     return s << std::string(d);
// //   }

// };


class RoiDescriptor {

  public:

    // iterator 
    typedef std::vector<const RoiDescriptor*>::const_iterator roi_iterator;
    /// convenient
    static constexpr bool FULLSCAN = true;
    static constexpr bool ROI = false;

   /**
   * @brief constructor
   * @param eta      eta of RoI
   * @param etaMinus eta at rear  of RoI
   * @param etaPlus  eta at front of RoI
   * @param phi      phi of RoI
   * @param phiMinus minimum phi of RoI
   * @param phiPlus  maximum phi of RoI
   * @param zed      zed of RoI
   * @param zedMinus zed at rear  of RoI
   * @param zedPlus  zed at front of RoI
   */
    RoiDescriptor(double eta,   double etaMinus,   double etaPlus, 
      double phi,   double phiMinus,   double phiPlus, 
      double zed=0, double zedMinus=-s_zedWidthDefault, double zedPlus=s_zedWidthDefault   );
      // zedminus - s_zedWidthDefault = 225 //from ROIDescriptor 

    /*
    *  need an explicit class copy constructor
    */
    RoiDescriptor( const RoiDescriptor& roi );
    RoiDescriptor& operator=( const RoiDescriptor& r );  

    
    // Destructor
    ~RoiDescriptor();


    // Methods to retrieve data members

    double phi() const  { return m_phi; }
    double eta() const  { return m_eta; }
    double zed() const { return m_zed; }

    /// these quantities probably don't need to be used any more
    /// - they are implemented here only because we had them in 
    ///   the original legacy interface

    double zedPlus()  const   { return m_zedPlus; } //!< z at the most forward end of the RoI
    double zedMinus() const   { return m_zedMinus; } //!< z at the most backward end of the RoI

    double etaPlus()  const { return m_etaPlus; }   //!< gets eta at zedPlus
    double etaMinus() const { return m_etaMinus; }   //!< gets eta at zMinus

    double phiPlus()  const { return m_phiPlus; }    //!< gets phiPlus
    double phiMinus() const { return m_phiMinus; }   //!< gets phiMinus


    /// versioning 
    int version() const { return m_version; }
    void version(int v)                  { m_version = v; }


    /// output
    // virtual operator std::string() const ;


    /// is this a full scan RoI?
    bool  isFullscan() const { return m_fullscan; }
  
    /// SuperRoI compatability methods

    /// am I a SuperRoi?
    bool composite() const { return m_composite; }
    void setComposite(bool b=true)          { m_composite=b; }

    /// always manage constituents ???
    bool manageConstituents() const { return m_manageConstituents; }
    void manageConstituents(bool b) { m_manageConstituents=b; }

    /// number of constituents
    unsigned size() const { return m_roiDescriptors.size(); }

    /// find an RoiDescriptor constituent
    const RoiDescriptor* at(int i) const { return m_roiDescriptors.at(i); }

    /// clear the vector
    void clear()  { m_roiDescriptors.clear(); }  // setComposite(false); }

    /// reserve elements in vector
    void reserve(size_t s) { m_roiDescriptors.reserve(s); }

    /// add a RoiDescriptor
    void push_back(const RoiDescriptor* roi) { m_roiDescriptors.push_back(roi); setComposite(true); }

    /// iterators
    roi_iterator  begin() const  { return m_roiDescriptors.begin(); }
    roi_iterator  end()   const  { return m_roiDescriptors.end(); }

    /// return the gradients
    double dzdrMinus() const { return m_dzdrMinus; }       //!<  dz/dr at the rear of the RoI
    double dzdrPlus()  const  { return m_dzdrPlus; }        //!<  dz/dr at the front of the RoI

    double drdzMinus() const { return m_drdzMinus; }       //!<  dr/dz at the rear of the RoI
    double drdzPlus()  const { return m_drdzPlus; }        //!<  dr/dz at the front of the RoI

    /// methods to calculate z position at the RoI boundary 
    /// at a given radius
    double zedMin(double r) const ;
    double zedMax(double r) const ;

    double zedOuterPlus()  const   { return m_zedOuterPlus; } //!< z at the most forward end of the RoI
    double zedOuterMinus() const  { return m_zedOuterMinus; } //!< z at the most backward end of the RoI

    double rhoMin(double z) const  ;
    double rhoMax(double z) const  ;

    static double zedWidthDefault() { return s_zedWidthDefault; }

    /// set default z-width (but only before any RoiDescriptor has been created)
    static void zedWidthDefault( double d );

    //fromn trig 
    unsigned int roiId() const   { return m_roiId; }
    unsigned int l1Id() const { return m_l1Id; }
    unsigned int roiWord() const  { return m_roiWord; }


  private:

    /// default parameters - there may be better ways, but this will do
    static std::atomic<double> s_zedWidthDefault;
    /// to ensure default width is only set once at job startup
    static std::atomic<bool> s_firstInstanceCreated;

    float m_phi;                 //!< phi of RoI center
    float m_eta;                 //!< eta of RoI center
    float m_zed;                 //!< zed of RoI center

    float m_phiMinus;            //!< most negative RoI in azimuthal
    float m_phiPlus;             //!< most positive RoI in azimuthal 
    float m_etaMinus;             //!< eta of RoI at zedMinus
    float m_etaPlus;              //!< eta of RoI at zedPlus
    float m_zedMinus;             //!< z position at most negative position along the beamline
    float m_zedPlus;              //!< z position at most positive position along the beamline

    float m_dzdrMinus;       //!<  dz/dr at the rear of the RoI
    float m_dzdrPlus;        //!<  dz/dr at the front of the RoI

    float m_drdzMinus;    //!<  dr/dz at the rear of the RoI
    float m_drdzPlus;     //!<  dr/dz at the front of the RoI

    float m_zedOuterMinus;  //!< z at rear of RoI at the outer radius ( = 1100 mm) 
    float m_zedOuterPlus;   //!< z at front of RoI at the outer radius ( = 1100 mm) 

    bool m_fullscan;             //!< flag this as a full detector RoI
    bool m_composite;            //!< flag this as a composite RoI
    bool m_manageConstituents;   //!< flag to determine whether consituents should be managed

    int  m_version;              //!< transient version identifier

    std::vector<const RoiDescriptor*>  m_roiDescriptors;  //!< roi constituents
    
    //from trig 
    unsigned int m_l1Id;          //!< lvl1 event number
    unsigned int m_roiId;         //!< RoI number
    unsigned int m_roiWord;       //!< lvl1 RoI word from which this RoI was initially constructed


//   std::string str( const RoiDescriptor& d );                           //<! printing helper
//   std::ostream& operator<<( std::ostream& m, const RoiDescriptor& d ); //<! printing helper (wraps above)

};






} //end of namesspace acts 

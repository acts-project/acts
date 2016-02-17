# A Tracking Software (ATS) Project

This library is based on the track reconstruction software developed by the 
[ATLAS Collaboration](http://cern.ch/atlas).

The main philosophy is to provide high-level track reconstruction modules that 
can be specified for detector technologies by simple extension.

* Event Data Model (EDM)
* Geometry 
* track reconstruction tools

The library is attempted to build against Gaudi and Gaudi-Athena, while
additional external dependencies are kept at a minimum.

# Gaudi build with CMake

source init.sh
(make clean) 
make -j 8
make install (to install library)


# Athena build with cmt

This is only possible with gcc >= 4.9

e.g. asetup,devval,rel0,opt,here,runtime
setupWorkArea
cd WorkArea/cmt
cmt br cmt config; cmt br gmake



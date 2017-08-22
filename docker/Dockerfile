FROM ubuntu:16.04

MAINTAINER Christian Gumpert <christian.gumpert@cern.ch>

SHELL ["/bin/bash","-c"]

# pre-requisites
RUN apt-get -y update && \
    apt-get -y --no-install-recommends install software-properties-common && \
    add-apt-repository -y ppa:ubuntu-toolchain-r/test && \
    apt-get -y update && \
    apt-get -y upgrade && \
    apt-get -y --no-install-recommends install \
      ca-certificates \
      clang-3.8 \
      cmake \
      doxygen \
      g++-6 \
      gcc-6 \
      git \
      graphviz \
      make \
      sudo \
      wget && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-6 20 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-6 20 && \
    update-alternatives --install /usr/bin/clang clang /usr/bin/clang-3.8 20 && \
    update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-3.8 20

# Eigen
RUN export EIGEN_VERSION="3.2.10" && \
    wget http://bitbucket.org/eigen/eigen/get/${EIGEN_VERSION}.tar.gz && \
    tar xf ${EIGEN_VERSION}.tar.gz && \
    mkdir -p /opt/ && \
    mv eigen-eigen-*/ /opt/eigen/ && \
    rm ${EIGEN_VERSION}.tar.gz

# boost
RUN export BOOST_VERSION="1.62.0" && export BOOST_VERSION_=${BOOST_VERSION//./_} && \
    wget -O boost-${BOOST_VERSION_}.tar.gz http://downloads.sourceforge.net/project/boost/boost/${BOOST_VERSION}/boost_${BOOST_VERSION_}.tar.gz && \
    tar xf boost-${BOOST_VERSION_}.tar.gz && \
    cd boost_${BOOST_VERSION_}/ && \
    ./bootstrap.sh --prefix=/opt/boost/ --with-libraries=program_options,test && \
    ./b2 install && \
    cd .. && rm boost-${BOOST_VERSION_}.tar.gz && rm -r boost_${BOOST_VERSION_}

# gosu
RUN export GOSU_VERSION="1.10" && \
    rm -rf /var/lib/apt/lists/* && \
    dpkgArch="$(dpkg --print-architecture | awk -F- '{ print $NF }')" && \
    wget -O /usr/local/bin/gosu "https://github.com/tianon/gosu/releases/download/${GOSU_VERSION}/gosu-$dpkgArch" && \
    wget -O /usr/local/bin/gosu.asc "https://github.com/tianon/gosu/releases/download/${GOSU_VERSION}/gosu-$dpkgArch.asc" && \
    export GNUPGHOME="$(mktemp -d)" && \
    gpg --keyserver ha.pool.sks-keyservers.net --recv-keys B42F6819007F00F88E364FD4036A9C25BF357DD4 && \
    gpg --batch --verify /usr/local/bin/gosu.asc /usr/local/bin/gosu && \
    rm -r "$GNUPGHOME" /usr/local/bin/gosu.asc && \
    chmod +x /usr/local/bin/gosu && \
    gosu nobody true
    
# ROOT
RUN export ROOT_VERSION="6.08.00" && \
    apt-get -y update && \
    apt-get -y upgrade && \
    apt-get -y --no-install-recommends install \
      libfftw3-dev \
      libgsl-dev \
      libx11-dev \
      libxext-dev \
      libxft-dev \
      libxml2-dev \
      libxpm-dev \
      python-dev \
      xlibmesa-glu-dev && \
    wget https://root.cern.ch/download/root_v${ROOT_VERSION}.source.tar.gz && \
    tar xf root_v${ROOT_VERSION}.source.tar.gz && \
    mkdir build && cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=/opt/root/ \
    	  -Dcxx14=ON \
	  -Dfftw3=ON \
	  -Dgdml=ON \
	  -Dminuit2=ON \
	  -Dopengl=ON \
	  -Droofit=ON \
	  -Dxml=ON \
	  ../root-${ROOT_VERSION}/ && \
    cmake --build . --target install -- -j 4 && \
    cd .. && rm -rf build/ && rm -rf root-${ROOT_VERSION}/ && rm root_v${ROOT_VERSION}.source.tar.gz

# CLHEP
RUN export CLHEP_VERSION="2.3.3.2" && \
    wget http://proj-clhep.web.cern.ch/proj-clhep/DISTRIBUTION/tarFiles/clhep-${CLHEP_VERSION}.tgz && \
    tar xf clhep-${CLHEP_VERSION}.tgz && \
    mkdir -p build && cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=/opt/clhep/ \
    	  -DCLHEP_BUILD_CXXSTD="-std=c++14" \
	  ../${CLHEP_VERSION}/CLHEP/ && \
    make -j 4 && \
    make install && \
    cd .. && rm -rf build/ && rm -rf ${CLHEP_VERSION}/ && rm clhep-${CLHEP_VERSION}.tgz

# GEANT4
RUN export GEANT4_VERSION="4.10.02.p02" && \
    apt-get update && \
    apt-get -y upgrade && \
    apt-get -y --no-install-recommends install libxerces-c-dev && \    
    wget http://geant4.cern.ch/support/source/geant${GEANT4_VERSION}.tar.gz && \
    tar xf geant${GEANT4_VERSION}.tar.gz && \
    mkdir -p build && cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=/opt/geant4/ \
    	  -DGEANT4_BUILD_CXXSTD=14 \
	  -DGEANT4_USE_GDML=ON \
	  -DCLHEP_ROOT_DIR=/opt/clhep/ \
	  ../geant${GEANT4_VERSION}/ && \
    make -j 4 && \
    make install && \
    cd .. && rm -rf build/ && rm -rf geant${GEANT4_VERSION}/ && rm geant${GEANT4_VERSION}.tar.gz

# DD4HEP
RUN export DD4HEP_VERSION="01-02" && \
    wget http://github.com/AIDASoft/DD4hep/archive/v${DD4HEP_VERSION}.tar.gz && \
    tar xf v${DD4HEP_VERSION}.tar.gz && \
    mkdir -p build && cd build && \
    export BOOST_ROOT=/opt/boost/ && \
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/root/lib:/opt/clhep/lib && \
    export PATH=$PATH:/opt/root/bin && \
    cmake -DCMAKE_INSTALL_PREFIX=/opt/dd4hep/ \
    	  -DBoost_NO_BOOST_CMAKE=ON \
	  -DDD4HEP_USE_GEANT4=ON \
	  -DDD4HEP_USE_CXX14=ON \
	  -DCLHEP_ROOT_DIR=/opt/clhep \
	  -DGeant4_DIR=/opt/geant4/lib/Geant4-10.2.2/ \
	  ../DD4hep-${DD4HEP_VERSION}/ && \
    make -j 4 && \
    make install && \
    cd .. && rm -rf build/ && rm -rf DD4hep-${DD4HEP_VERSION}/ && rm v${DD4HEP_VERSION}.tar.gz
    
ENV PATH="$PATH:/opt/root/bin:/opt/clhep/bin:/opt/geant4/bin:/opt/dd4hep/bin"

# set working directory
WORKDIR /workdir

# set up start environment
COPY acts-build acts-start /usr/local/bin/

# start a bash shell by default
ENTRYPOINT ["/usr/local/bin/acts-start"]
CMD ["/bin/bash"]

# setup appropriate LCG 100 release via cvmfs

if test -e /etc/centos-release && grep 'CentOS Linux release 7' /etc/centos-release; then
  lcg_os=centos7
elif test -e /etc/centos-release && grep 'CentOS Stream release 8' /etc/centos-release; then
  lcg_os=centos8
else
  echo "Unsupported system" 1>&2
  return
fi

lcg_release=LCG_100
lcg_compiler=gcc10-opt
lcg_platform=x86_64-${lcg_os}-${lcg_compiler}
lcg_view=/cvmfs/sft.cern.ch/lcg/views/${lcg_release}/${lcg_platform}

source ${lcg_view}/setup.sh
# extra variables required to build acts
export DD4hep_DIR=${lcg_view}

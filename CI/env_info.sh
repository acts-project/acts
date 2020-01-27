if [ "$(cat /etc/redhat-release | grep '^Scientific Linux .*6.*')" ]; then
  export ACTS_OS=slc6
  export ACTS_COMPILER=gcc8-opt
elif [ "$(cat /etc/centos-release | grep 'CentOS Linux release 7')" ]; then
  export ACTS_OS=centos7
  export ACTS_COMPILER=gcc8-opt
elif [ "$(cat /etc/debian_version | grep 'buster.*sid')" ]; then
  export ACTS_OS=ubuntu1804
  if [[ "${1}" =~ "lcg94" ]]; then
    export ACTS_COMPILER=gcc7-opt
  else
    export ACTS_COMPILER=gcc8-opt
  fi
else
  echo "Unknown OS" 1>&2
fi

echo "$ACTS_OS $ACTS_COMPILER"

if [ "$(echo $CI_RUNNER_TAGS | grep 'acts')" ]; then
  echo "Running on Acts CI runner: use 8 cores"
  export ACTS_NCPUS=8
else
  echo "NOT running on Acts CI runner: use 2 cores"
  export ACTS_NCPUS=2
fi


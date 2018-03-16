if [ "$(cat /etc/redhat-release | grep '^Scientific Linux .*6.*')" ]; then
  export ACTS_OS=slc6
  export ACTS_COMPILER=gcc62-opt
elif [ "$(cat /etc/centos-release | grep 'CentOS Linux release 7')" ]; then
  export ACTS_OS=centos7
  export ACTS_COMPILER=gcc7-opt
else
  echo "Unknown OS" 1>&2
  exit 1
fi

echo "$ACTS_OS $ACTS_COMPILER"

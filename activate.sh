. ~/cern/source/spack/share/spack/setup-env.sh
spack env activate ~/cern/source/acts/ci-dependencies
. ~/cern/source/acts/ci-dependencies/.python-env/bin/activate
. ~/cern/build/acts/acts/dev2/this_acts_withdeps.sh

#TODO include into `this_acts_withdeps.sh`
export PYTHONPATH=~/cern/source/acts/acts/dev2/Examples/Scripts/Python:$PYTHONPATH
export PYTHONPATH=~/cern/source/acts/acts/dev2/CI/physmon:$PYTHONPATH

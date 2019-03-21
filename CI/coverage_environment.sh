SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
dir=$(mktemp -d)

target_command=$1
echo $target_command

function finish {
  rm -rf $dir
}
trap finish EXIT

python=/cvmfs/sft.cern.ch/lcg/views/LCG_94python3/x86_64-slc6-gcc8-opt/bin/python

$python -m venv $dir

source $dir/bin/activate

which python
which pip

echo $SCRIPTPATH
pip install -r $SCRIPTPATH/requirements.txt

exec bash -c "$target_command"

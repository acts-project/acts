SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}"  )" &> /dev/null && pwd  )
export PYTHONPATH+=:$(realpath "${SCRIPT_DIR}/../../Examples/Scripts/Python/")
export PYTHONPATH+=":${SCRIPT_DIR}"

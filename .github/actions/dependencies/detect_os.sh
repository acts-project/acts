if [ $(uname) == "Linux" ]; then
  os_name=$(cat /etc/os-release | grep -e "^PRETTY_NAME=" | sed 's/PRETTY_NAME="\(.*\)"/\1/g')
  if [[ $os_name == *"Ubuntu"* ]]; then
    os="ubuntu"
  elif [[ $os_name == *"AlmaLinux"* ]]; then
    os="almalinux"
  fi
elif [ $(uname) == "Darwin" ]; then
  os_name="$(sw_vers -productName) $(sw_vers -productVersion)"
  os="macos"
else
  echo "Only Ubuntu, AlmaLinux and macOS are supported. Exiting."
  exit 1
fi

export os
export os_name

current_tag=$(git describe --exact-match HEAD 2>&1)
if [ "$?" -ne 0 ]; then
  echo "Not on tagged commit"
  # use previous tag on branch
  export prev_tag=$(git describe --tags --abbrev=0)
else
  echo "On tagged commit $current_tag"
  # use tag before this tag
  export prev_tag=$(git describe --tags --abbrev=0 HEAD~)
fi

echo "Previous tag to build first $prev_tag"

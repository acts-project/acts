if [ "${BASH_SOURCE[0]}" -ef "$0" ]; then
    echo "This script MUST be sourced, not run!"
    exit 1
fi

CHANGESET_FILE=${GITHUB_WORKSPACE}/.changeset.txt

if grep "^Core/" ${CHANGESET_FILE}; then
  echo "touches_core=true" >> "$GITHUB_OUTPUT"
else
  echo "touches_core=false" >> "$GITHUB_OUTPUT"
fi

if grep "^Examples/" ${CHANGESET_FILE}; then
  echo "touches_examples=true" >> "$GITHUB_OUTPUT"
else
  echo "touches_examples=false" >> "$GITHUB_OUTPUT"
fi

if grep "^CI/" ${CHANGESET_FILE}; then
  echo "touches_infra=true" >> "$GITHUB_OUTPUT"
else
  echo "touches_infra=false" >> "$GITHUB_OUTPUT"
fi

if grep "^Detray/" ${CHANGESET_FILE}; then
  echo "touches_detray=true" >> "$GITHUB_OUTPUT"
else
  echo "touches_detray=false" >> "$GITHUB_OUTPUT"
fi

CHANGESET_FILE=$1
OUTPUT_FILE=$2

check_pattern () {
  if grep $1 ${CHANGESET_FILE}; then
    echo "TOUCHES_$2=true" >> "$OUTPUT_FILE"
  else
    echo "TOUCHES_$2=false" >> "$OUTPUT_FILE"
  fi
}

check_pattern "^Core/" "CORE"
check_pattern "^Examples/" "EXAMPLES"
check_pattern "^Detray/" "DETRAY"
check_pattern "^Traccc/" "TRACCC"

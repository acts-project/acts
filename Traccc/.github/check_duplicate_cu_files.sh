#!/bin/bash

! [[ $(find device tests -name "*\.cu" | xargs -I {} basename {} | sort | uniq -d) ]]
exit $?

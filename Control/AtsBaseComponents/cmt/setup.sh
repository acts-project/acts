# echo "setup AtsBaseComponents AtsBaseComponents-00-06-13-01 in /afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/a-tracking-sw/Control"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/atlas.cern.ch/repo/sw/software/x86_64-slc6-gcc48-opt/20.7.1/CMT/v1r25p20140131; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtAtsBaseComponentstempfile=`${CMTROOT}/${CMTBIN}/cmt.exe -quiet build temporary_name`
if test ! $? = 0 ; then cmtAtsBaseComponentstempfile=/tmp/cmt.$$; fi
${CMTROOT}/${CMTBIN}/cmt.exe setup -sh -pack=AtsBaseComponents -version=AtsBaseComponents-00-06-13-01 -path=/afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/a-tracking-sw/Control  -no_cleanup $* >${cmtAtsBaseComponentstempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/${CMTBIN}/cmt.exe setup -sh -pack=AtsBaseComponents -version=AtsBaseComponents-00-06-13-01 -path=/afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/a-tracking-sw/Control  -no_cleanup $* >${cmtAtsBaseComponentstempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtAtsBaseComponentstempfile}
  unset cmtAtsBaseComponentstempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtAtsBaseComponentstempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtAtsBaseComponentstempfile}
unset cmtAtsBaseComponentstempfile
return $cmtsetupstatus


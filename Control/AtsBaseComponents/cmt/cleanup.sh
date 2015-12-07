# echo "cleanup AtsBaseComponents AtsBaseComponents-00-06-13-01 in /afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/a-tracking-sw/Control"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/atlas.cern.ch/repo/sw/software/x86_64-slc6-gcc48-opt/20.7.1/CMT/v1r25p20140131; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtAtsBaseComponentstempfile=`${CMTROOT}/${CMTBIN}/cmt.exe -quiet build temporary_name`
if test ! $? = 0 ; then cmtAtsBaseComponentstempfile=/tmp/cmt.$$; fi
${CMTROOT}/${CMTBIN}/cmt.exe cleanup -sh -pack=AtsBaseComponents -version=AtsBaseComponents-00-06-13-01 -path=/afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/a-tracking-sw/Control  $* >${cmtAtsBaseComponentstempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/${CMTBIN}/cmt.exe cleanup -sh -pack=AtsBaseComponents -version=AtsBaseComponents-00-06-13-01 -path=/afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/a-tracking-sw/Control  $* >${cmtAtsBaseComponentstempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtAtsBaseComponentstempfile}
  unset cmtAtsBaseComponentstempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtAtsBaseComponentstempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtAtsBaseComponentstempfile}
unset cmtAtsBaseComponentstempfile
return $cmtcleanupstatus


# echo "setup EventPrimitives EventPrimitives-00-00-45 in /afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/a-tracking-sw/Event"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/atlas.cern.ch/repo/sw/software/x86_64-slc6-gcc48-opt/20.7.1/CMT/v1r25p20140131; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtEventPrimitivestempfile=`${CMTROOT}/${CMTBIN}/cmt.exe -quiet build temporary_name`
if test ! $? = 0 ; then cmtEventPrimitivestempfile=/tmp/cmt.$$; fi
${CMTROOT}/${CMTBIN}/cmt.exe setup -sh -pack=EventPrimitives -version=EventPrimitives-00-00-45 -path=/afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/a-tracking-sw/Event  -no_cleanup $* >${cmtEventPrimitivestempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/${CMTBIN}/cmt.exe setup -sh -pack=EventPrimitives -version=EventPrimitives-00-00-45 -path=/afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/a-tracking-sw/Event  -no_cleanup $* >${cmtEventPrimitivestempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtEventPrimitivestempfile}
  unset cmtEventPrimitivestempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtEventPrimitivestempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtEventPrimitivestempfile}
unset cmtEventPrimitivestempfile
return $cmtsetupstatus


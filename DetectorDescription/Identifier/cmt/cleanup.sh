# echo "cleanup Identifier Identifier-00-09-32 in /afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/a-tracking-sw/DetectorDescription"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/atlas.cern.ch/repo/sw/software/x86_64-slc6-gcc48-opt/20.7.1/CMT/v1r25p20140131; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtIdentifiertempfile=`${CMTROOT}/${CMTBIN}/cmt.exe -quiet build temporary_name`
if test ! $? = 0 ; then cmtIdentifiertempfile=/tmp/cmt.$$; fi
${CMTROOT}/${CMTBIN}/cmt.exe cleanup -sh -pack=Identifier -version=Identifier-00-09-32 -path=/afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/a-tracking-sw/DetectorDescription  $* >${cmtIdentifiertempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/${CMTBIN}/cmt.exe cleanup -sh -pack=Identifier -version=Identifier-00-09-32 -path=/afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/a-tracking-sw/DetectorDescription  $* >${cmtIdentifiertempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtIdentifiertempfile}
  unset cmtIdentifiertempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtIdentifiertempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtIdentifiertempfile}
unset cmtIdentifiertempfile
return $cmtcleanupstatus


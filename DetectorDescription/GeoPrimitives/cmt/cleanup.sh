# echo "cleanup GeoPrimitives GeoPrimitives-00-00-33-05 in /afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw/DetectorDescription"

if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/cern.ch/sw/contrib/CMT/v1r25p20140131; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtGeoPrimitivestempfile=`${CMTROOT}/${CMTBIN}/cmt.exe -quiet build temporary_name`
if test ! $? = 0 ; then cmtGeoPrimitivestempfile=/tmp/cmt.$$; fi
${CMTROOT}/${CMTBIN}/cmt.exe cleanup -sh -pack=GeoPrimitives -version=GeoPrimitives-00-00-33-05 -path=/afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw/DetectorDescription  $* >${cmtGeoPrimitivestempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/${CMTBIN}/cmt.exe cleanup -sh -pack=GeoPrimitives -version=GeoPrimitives-00-00-33-05 -path=/afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw/DetectorDescription  $* >${cmtGeoPrimitivestempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtGeoPrimitivestempfile}
  unset cmtGeoPrimitivestempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtGeoPrimitivestempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtGeoPrimitivestempfile}
unset cmtGeoPrimitivestempfile
return $cmtcleanupstatus


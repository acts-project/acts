# echo "setup GeoPrimitives GeoPrimitives-00-00-33-05 in /afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw/DetectorDescription"

if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/cern.ch/sw/contrib/CMT/v1r25p20140131; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtGeoPrimitivestempfile=`${CMTROOT}/${CMTBIN}/cmt.exe -quiet build temporary_name`
if test ! $? = 0 ; then cmtGeoPrimitivestempfile=/tmp/cmt.$$; fi
${CMTROOT}/${CMTBIN}/cmt.exe setup -sh -pack=GeoPrimitives -version=GeoPrimitives-00-00-33-05 -path=/afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw/DetectorDescription  -no_cleanup $* >${cmtGeoPrimitivestempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/${CMTBIN}/cmt.exe setup -sh -pack=GeoPrimitives -version=GeoPrimitives-00-00-33-05 -path=/afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw/DetectorDescription  -no_cleanup $* >${cmtGeoPrimitivestempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtGeoPrimitivestempfile}
  unset cmtGeoPrimitivestempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtGeoPrimitivestempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtGeoPrimitivestempfile}
unset cmtGeoPrimitivestempfile
return $cmtsetupstatus


# echo "setup Identifier Identifier-00-09-32 in /afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw/DetectorDescription"

if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/cern.ch/sw/contrib/CMT/v1r25p20140131; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtIdentifiertempfile=`${CMTROOT}/${CMTBIN}/cmt.exe -quiet build temporary_name`
if test ! $? = 0 ; then cmtIdentifiertempfile=/tmp/cmt.$$; fi
${CMTROOT}/${CMTBIN}/cmt.exe setup -sh -pack=Identifier -version=Identifier-00-09-32 -path=/afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw/DetectorDescription  -no_cleanup $* >${cmtIdentifiertempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/${CMTBIN}/cmt.exe setup -sh -pack=Identifier -version=Identifier-00-09-32 -path=/afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw/DetectorDescription  -no_cleanup $* >${cmtIdentifiertempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtIdentifiertempfile}
  unset cmtIdentifiertempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtIdentifiertempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtIdentifiertempfile}
unset cmtIdentifiertempfile
return $cmtsetupstatus


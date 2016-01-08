# echo "cleanup Identifier Identifier-00-09-32 in /afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw/DetectorDescription"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /afs/cern.ch/sw/contrib/CMT/v1r25p20140131
endif
source ${CMTROOT}/mgr/setup.csh
set cmtIdentifiertempfile=`${CMTROOT}/${CMTBIN}/cmt.exe -quiet build temporary_name`
if $status != 0 then
  set cmtIdentifiertempfile=/tmp/cmt.$$
endif
${CMTROOT}/${CMTBIN}/cmt.exe cleanup -csh -pack=Identifier -version=Identifier-00-09-32 -path=/afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw/DetectorDescription  $* >${cmtIdentifiertempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/${CMTBIN}/cmt.exe cleanup -csh -pack=Identifier -version=Identifier-00-09-32 -path=/afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw/DetectorDescription  $* >${cmtIdentifiertempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtIdentifiertempfile}
  unset cmtIdentifiertempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtIdentifiertempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtIdentifiertempfile}
unset cmtIdentifiertempfile
exit $cmtcleanupstatus


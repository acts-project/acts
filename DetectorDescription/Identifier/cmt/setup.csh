# echo "setup Identifier Identifier-00-09-32 in /afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/a-tracking-sw/DetectorDescription"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/atlas.cern.ch/repo/sw/software/x86_64-slc6-gcc48-opt/20.7.1/CMT/v1r25p20140131
endif
source ${CMTROOT}/mgr/setup.csh
set cmtIdentifiertempfile=`${CMTROOT}/${CMTBIN}/cmt.exe -quiet build temporary_name`
if $status != 0 then
  set cmtIdentifiertempfile=/tmp/cmt.$$
endif
${CMTROOT}/${CMTBIN}/cmt.exe setup -csh -pack=Identifier -version=Identifier-00-09-32 -path=/afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/a-tracking-sw/DetectorDescription  -no_cleanup $* >${cmtIdentifiertempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/${CMTBIN}/cmt.exe setup -csh -pack=Identifier -version=Identifier-00-09-32 -path=/afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/a-tracking-sw/DetectorDescription  -no_cleanup $* >${cmtIdentifiertempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtIdentifiertempfile}
  unset cmtIdentifiertempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtIdentifiertempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtIdentifiertempfile}
unset cmtIdentifiertempfile
exit $cmtsetupstatus


# echo "cleanup AtsBaseComponents AtsBaseComponents-00-06-13-01 in /afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/a-tracking-sw/Control"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/atlas.cern.ch/repo/sw/software/x86_64-slc6-gcc48-opt/20.7.1/CMT/v1r25p20140131
endif
source ${CMTROOT}/mgr/setup.csh
set cmtAtsBaseComponentstempfile=`${CMTROOT}/${CMTBIN}/cmt.exe -quiet build temporary_name`
if $status != 0 then
  set cmtAtsBaseComponentstempfile=/tmp/cmt.$$
endif
${CMTROOT}/${CMTBIN}/cmt.exe cleanup -csh -pack=AtsBaseComponents -version=AtsBaseComponents-00-06-13-01 -path=/afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/a-tracking-sw/Control  $* >${cmtAtsBaseComponentstempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/${CMTBIN}/cmt.exe cleanup -csh -pack=AtsBaseComponents -version=AtsBaseComponents-00-06-13-01 -path=/afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/a-tracking-sw/Control  $* >${cmtAtsBaseComponentstempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtAtsBaseComponentstempfile}
  unset cmtAtsBaseComponentstempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtAtsBaseComponentstempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtAtsBaseComponentstempfile}
unset cmtAtsBaseComponentstempfile
exit $cmtcleanupstatus


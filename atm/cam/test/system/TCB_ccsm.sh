#!/bin/sh 
#

if [ $# -ne 2 ]; then
    echo "TCB_ccsm.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TCB_ccsm.$1.$2

if [ -f ${CAM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CAM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TCB_ccsm.sh: ccsm configure and build test has already passed; results are in "
	echo "        ${CAM_TESTDIR}/${test_name}" 
        exit 0
    else
	read fail_msg < ${CAM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TCB_ccsm.sh: ccsm configure and build test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CAM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TCB_ccsm.sh: this ccsm configure and build test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CAM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CAM_TESTDIR}/${test_name} ${CAM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi

blddir=${CAM_TESTDIR}/${test_name}
if [ -d ${blddir} ]; then
    rm -rf ${blddir}
fi
mkdir -p ${blddir} 
if [ $? -ne 0 ]; then
    echo "TCB_ccsm.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${blddir}

echo "TCB_ccsm.sh: building ccsm executable; output in ${CAM_TESTDIR}/${test_name}/test.log" 
if [ -d ${CAM_TESTDIR}/case.$1.$2 ]; then
    rm -rf ${CAM_TESTDIR}/case.$1.$2
fi

${CAM_ROOT}/scripts/create_newcase -case ${CAM_TESTDIR}/case.$1.$2 -res $1 -compset $2 -mach ${CCSM_MACH} > test.log 2>&1
rc=$?
if [ $rc -eq 0 ]; then
    echo "TCB_ccsm.sh: create_newcase was successful" 
else
    echo "TCB_ccsm.sh: create_newcase failed, error from create_newcase= $rc" 
    echo "TCB_ccsm.sh: see ${CAM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

cd ${CAM_TESTDIR}/case.$1.$2
./xmlchange -file env_build.xml -id EXEROOT -val ${CAM_TESTDIR}/case.$1.$2 -silent

cp ${CAM_SCRIPTDIR}/nl_files/user_nl_cam .

./configure -case >> ${CAM_TESTDIR}/${test_name}/test.log 2>&1
rc=$?
if [ $rc -eq 0 ]; then
    echo "TCB_ccsm.sh: ccsm configure was successful" 
else
    echo "TCB_ccsm.sh: ccsm configure failed, error from configure= $rc" 
    echo "TCB_ccsm.sh: see ${CAM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi

buildscript=`ls *.build`
./$buildscript >> ${CAM_TESTDIR}/${test_name}/test.log 2>&1
rc=$?
if [ $rc -eq 0 ]; then
    echo "TCB_ccsm.sh: ccsm build was successful" 
else
    echo "TCB_ccsm.sh: ccsm build failed, error from build= $rc" 
    echo "TCB_ccsm.sh: see ${CAM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

cd Buildconf
cat cam.buildnml.csh | sed '
/cam_inparm/ a\
 nhtfrq=-24 \
 mfilt=1
' > cam.buildnml.csh.tmp
chmod a+x cam.buildnml.csh.tmp
mv cam.buildnml.csh.tmp cam.buildnml.csh

cd ${blddir}
echo "TCB_ccsm.sh: ccsm configure and build test passed"
echo "PASS" > TestStatus
if [ $CAM_RETAIN_FILES != "TRUE" ]; then
    echo "TCB_ccsm.sh: removing some unneeded files to save disc space" 
    rm -rf ${CAM_TESTDIR}/case.$1.$2/atm
    rm -rf ${CAM_TESTDIR}/case.$1.$2/glc
    rm -rf ${CAM_TESTDIR}/case.$1.$2/ice
    rm -rf ${CAM_TESTDIR}/case.$1.$2/lnd
    rm -rf ${CAM_TESTDIR}/case.$1.$2/ocn
    rm -rf ${CAM_TESTDIR}/case.$1.$2/cpl
    rm -rf ${CAM_TESTDIR}/case.$1.$2/mct
    rm -rf ${CAM_TESTDIR}/case.$1.$2/pio
    rm -rf ${CAM_TESTDIR}/case.$1.$2/ccsm
    rm -rf ${CAM_TESTDIR}/case.$1.$2/csm_share
    rm -rf ${CAM_TESTDIR}/case.$1.$2/lib
fi

exit 0

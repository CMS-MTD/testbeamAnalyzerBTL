#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_89/x86_64-slc6-gcc62-opt/setup.sh
cd ${_CONDOR_SCRATCH_DIR}
#setenv USER `whoami`
export USER=$(whoami)
echo "User is"
echo $USER
tar -xf $6
rm $6

python analyze_base.py --bar ${1} --firstRun ${2} --lastRun ${3} --biasVoltage ${4} --timeAlgo ${5} --isCondor 

### Now that the run is over, there is one or more root files created
echo "List all produced .png files = "
ls $7/*.png
ls $7/*.txt
echo "List all files"
ls 
echo "*******************************************"
OUTDIR=root://cmseos.fnal.gov//store/user/$USER/testbeam_04-2019/
echo "xrdcp output for condor"

# copy .png files
for FILE in $7/*.png
do
  echo "xrdcp -f ${FILE} ${OUTDIR}/${FILE}"
  xrdcp -f ${FILE} ${OUTDIR}/${FILE} 2>&1
  XRDEXIT=$?
  if [[ $XRDEXIT -ne 0 ]]; then
    rm *.png
    echo "exit code $XRDEXIT, failure in xrdcp"
    exit $XRDEXIT
  fi
  rm ${FILE}
done

# copy .txt files
for FILE in $7/*.txt
do
  echo "xrdcp -f ${FILE} ${OUTDIR}/${FILE}"
  xrdcp -f ${FILE} ${OUTDIR}/${FILE} 2>&1
  XRDEXIT=$?
  if [[ $XRDEXIT -ne 0 ]]; then
    rm *.png
    echo "exit code $XRDEXIT, failure in xrdcp"
    exit $XRDEXIT
  fi
  rm ${FILE}
done

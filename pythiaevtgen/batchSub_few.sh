#!/bin/bash

[ "$USER" == "bpn7" ]     && WorkDir=/nfs/slac/g/atlas/u01/users/bnachman/SLAC_pythia/PURJ/src/

SubFileLoc=`pwd`/_batchSingleSub.sh
#rm $SubFileLoc
DateSuffix=`date +%Y%m%d_%Hh%Mmin`

echo '#!/bin/bash
echo CD to $1
echo CMD is $2

cd $1
source /nfs/slac/g/atlas/u01/users/bnachman/SLAC_pythia/PURJ/src/slac_specific/setup.sh
cmd=$4

echo MAKING TEMP DIR $2
JOBFILEDIR=$2
mkdir $JOBFILEDIR
REALOUT=$3
echo MADE TEMP DIR $JOBFILEDIR
echo WILL COPY TO $REALOUT

shift
shift
echo Calling $cmd $*
$cmd $*
echo ls $JOBFILEDIR
ls $JOBFILEDIR
cp -r $JOBFILEDIR/*.txt $REALOUT
echo COPYING to $REALOUT
rm -rf $JOBFILEDIR
' > $SubFileLoc
chmod u+x $SubFileLoc

#----------------
#1000 for less than 18.
#500 between 18 and 47.
#200 between 48 and 78.
#50 between 79 and 85.
#10 beyond that.
#Queue=long
Queue=atlas-t3
#Queue=long 
nevents=10
njobs=3
npu_min=$1
npu_max=$1
LogPrefix=`pwd`/logs/${DateSuffix}/${DateSuffix}_bsub_${mu}_
OutDirFinal=`pwd`/files/${DateSuffix}
mkdir -p `dirname $LogPrefix`
mkdir -p $OutDirFinal
echo
echo "Submitting $njobs jobs each with $nevents events to $Queue"
echo $LogPrefix
echo ${njobs},${npu_min}

for (( npu=$npu_min; npu<=$npu_max; npu++ )) ;  do
    for (( ii=1; ii<=$njobs; ii++ )) ;  do
	echo $ii
	OutDir=/scratch/${DateSuffix}_${ii}/
	bsub -q ${Queue} -R 'select[(!preempt&&rhel60&&cvmfs&&inet)]' -o $LogPrefix${ii}_${npu}_.log $SubFileLoc           \
            ${WorkDir} ${OutDir} ${OutDirFinal} ./dijets.exe  \
            -out ${OutDir}/Sample_mu_${npu}_nevents_${nevents}_job_${ii}.txt \
            -nev ${nevents} -npu ${npu} -jet-only
    done
done
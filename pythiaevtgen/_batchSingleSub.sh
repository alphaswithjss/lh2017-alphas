#!/bin/bash
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


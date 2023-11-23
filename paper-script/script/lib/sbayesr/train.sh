#!/bin/bash
#
# assumed to be called from ./script/nonaddpgs/boosting/pipeline/p1-2_boosting.sh
#

set -eux

program="sbayesr"

echo "Start $(readlink -f $0)"
hostname
if [ -v SGE_O_WORKDIR ]; then
	# assume you are in the root
	:
else
	# to root
	cd "$(dirname $0)/../../../"
fi
pwd

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh

if [ $# -ne 3 ]; then
	echo "Wrong argument number."
	exit 1
fi

dout="$1"
fma="$2"
fldm="$3"

echo "dout: $dout"
echo "fma: $fma"
echo "fldm: $fldm"

mkdir -p "$dout"

# TODO
flog="${dout}/${program}.$(date +'%y%m%d-%H%M%S').log"

date "+%Y-%m-%d %H:%M:%S"

if [ -e ${dout}out.snpRes ]; then
	echo "result exists"
	exit 0
fi

if [ ! -e ${fma} ]; then
	echo "fma does not exist"
	exit 0
fi

source ./script/lib/system/module.sh load gctb/2.02

gctb --sbayes R \
	--ldm ${fldm} \
	--gwas-summary ${fma} \
	--out ${dout}out \
	--p-value 0.4 \
	--chain-length 10000 --burn-in 2000 --out-freq 10 |&
	tee $flog

date "+%Y-%m-%d %H:%M:%S"

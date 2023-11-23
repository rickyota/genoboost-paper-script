#!/bin/bash

set -eux

program="ldpred"

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
source ./script/lib/system/conda_activate.sh $program

if [ $# -ne 6 ]; then
	echo "Wrong argument number."
	exit 1
fi

dout="$1"
fss="$2"
gtref="$3"
nss="$4"
fs="$5"
threads="$6"

mkdir -p "$dout"

# skip wgt
#if [ -e ${fout_rust}.wgt ];then
#    echo "wgt exists: ${fout_rust}.wgt"
#    exit 0
#fi

# use all var
maf=0.01

# > 1M/3k~300
ldr=500

niter=100
nburnin=10

resultcoord="${dout}/coord"
resultld="${dout}/ld"
resultwgt="${dout}/wgt"

echo "$resultcoord"

# absolute path
resultcoord=$(readlink -f "$resultcoord")
resultld=$(readlink -f "$resultld")
resultwgt=$(readlink -f "$resultwgt")
fss=$(readlink -f "$fss")
gtref=$(readlink -f "$gtref")

source ./script/lib/system/python_threads.sh $threads

date "+%Y-%m-%d %H:%M:%S"

cd ./lib/ldpred/v1.0.11

if [ ! -e $resultcoord ]; then
	python -m ldpred \
		--debug coord \
		--out $resultcoord --gf $gtref --ssf $fss \
		--maf $maf --N $nss \
		--eff_type OR --ssf-format CUSTOM --ncol OBS_CT \
		--rs ID --A1 A1 --A2 A2 --pos POS --chr CHROM --pval P --eff OR
else
	echo "coord exist: $resultcoord"
fi

if [ ! -e "${resultwgt}_LDpred_p1.0000e-01.txt" ]; then
	python -m ldpred \
		--debug gibbs \
		--out $resultwgt --cf $resultcoord --ldf $resultld \
		--ldr $ldr --f $fs --N $nss --n-iter $niter --n-burn-in $nburnin
fi
date "+%Y-%m-%d %H:%M:%S"

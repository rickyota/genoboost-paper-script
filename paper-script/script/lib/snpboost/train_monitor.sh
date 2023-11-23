#!/bin/bash

set -eux

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh

if [ $# -ne 8 ]; then
	echo "Wrong argument number."
	exit 1
fi

phe="$1"
dout="$2"
fgenot="$3"
fphe="$4"
fsamples="$5"
sex="$6"
cvi="$7"
threads="$8"

mkdir -p "$dout"

# env
if [[ $HPCL_HOST = "pg" ]]; then
	set +ux
	source /bio/package/R/setup-4.0.3.sh
	set -ux
	source ./script/lib/system/module.sh load plink2
elif [[ $HPCL_HOST = "obcx" ]]; then
	source ./script/lib/system/module.sh load gcc/7.5.0
	PATH="/work/gg57/j29002/.local/bin/R-4.0.3/bin:$PATH"
	PATH="/work/gg57/j29002/.local/bin/plink-2.0:$PATH"
	PATH="/work/gg57/j29002/.local/bin/zstd-1.5.2:$PATH"
fi

# library path to pass
if [[ $HPCL_HOST = "pg" ]]; then
	libpath="/nfs/data06/ricky/local/lib/snpnet_R4.0.3"
elif [[ $HPCL_HOST = "obcx" ]]; then
	libpath="/work/gg57/j29002/.local/lib/R-4.0.3_snpnet2"
fi

date "+%Y-%m-%d %H:%M:%S"

if [ ! -f ${dout}/out.out.rds ]; then

	if [ ! -f ${dout}out/results/output_iter_1.RData ]; then
		echo "run from scratch"

		# for Rfast
		export USE_OPENMP=1
		# TEMP: otherwise warning keeps printing; OpenBLAS Warning : Detect OpenMP Loop and this application may hang. Please rebuild the library with USE_OPENMP=1 option.
		export OMP_NUM_THREADS=1
		Rscript ./lib/snpboost/v230912/src/run_cov_monitor.R \
			$libpath \
			${dout}/out $fgenot $fphe $phe \
			$fsamples \
			$sex $cvi \
			$threads

	else
		echo "ny"
	fi

fi

date "+%Y-%m-%d %H:%M:%S"

if [ -f ${dout}/out.out.rds ]; then
	# if out file exists
	# extract wgt
	Rscript ./lib/snpboost/v230912/src/run_loadrd.R $dout
else
	echo "better not to use."
	exit 1
fi

date "+%Y-%m-%d %H:%M:%S"

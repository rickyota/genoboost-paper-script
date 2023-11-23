#!/bin/bash
#
# assumed to be called from ./script/nonaddpgs/boosting/pipeline/p1-2_boosting.sh
#

set -eux

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

if [ $# -ne 7 ]; then
	echo "Wrong argument number."
	exit 1
fi

phe="$1"
dout="$2"
fgenot="$3"
fphe="$4"
fsamples="$5"
sex="$6"
threads="$7"

echo "dout: $dout"
echo "fgenot: $fgenot"
echo "fphe: $fphe"
echo "phe: $phe"
echo "sex: $sex"
echo "fsamples: $fsamples"
echo "threads: $threads"

mkdir -p "$dout"

# env
set +ux
if [[ $HPCL_HOST = "pg" ]]; then
	source /bio/package/R/setup-4.0.3.sh
	source ./script/lib/system/module.sh load plink2
elif [[ $HPCL_HOST = "obcx" ]]; then
	source ./script/lib/system/module.sh load gcc/7.5.0
	PATH="/work/gg57/j29002/.local/bin/R-4.0.3/bin:$PATH"
	PATH="/work/gg57/j29002/.local/bin/plink-2.0:$PATH"
	PATH="/work/gg57/j29002/.local/bin/zstd-1.5.2:$PATH"
fi
set -ux

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

		Rscript ./lib/snpnet/v210529/src/run_cov_snpnet2.R \
			$libpath \
			${dout}/out $fgenot $fphe $phe \
			$fsamples \
			$sex $threads
	else
		echo "ny"
		exit 1
	fi

fi

date "+%Y-%m-%d %H:%M:%S"

if [ -f ${dout}/out.out.rds ]; then
	# if out file exists
	# extract wgt
	Rscript ./lib/snpnet/v210529/src/run_loadrd.R $dout
else
	echo "ny"
	exit 1
fi

date "+%Y-%m-%d %H:%M:%S"

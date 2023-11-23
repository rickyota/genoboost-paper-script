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

if [ $# -ne 4 ]; then
	echo "Wrong argument number."
	exit 1
fi

dout="$1"
fss="$2"
gtref="$3"
nss="$4"

echo "dout: $dout"
echo "fss: $fss"
echo "gtref: $gtref"
echo "nss: $nss"

mkdir -p "$dout"

# TODO: should be above
flog="${dout}/snpnet.$(date +'%y%m%d-%H%M%S').log"
fout="${dout}/lassosum"

# env
set +ux
if [[ $HPCL_HOST = "pg" ]]; then
	source /bio/package/R/setup-4.0.3.sh
	#source /bio/package/R/setup-3.6.1.sh
	export PATH="/home/ricky/.local/lib/R4.0.3/lassosum/:$PATH"
elif [[ $HPCL_HOST = "obcx" ]]; then
	echo "ny"
	exit 1
	source ./script/lib/system/module.sh load gcc/7.5.0
	PATH="/work/gg57/j29002/.local/bin/R-4.0.3/bin:$PATH"
	#PATH="/work/gg57/j29002/.local/bin/plink-2.0:$PATH"
	#PATH="/work/gg57/j29002/.local/bin/zstd-1.5.2:$PATH"
fi
set -ux

## library path to pass
#if [[ $HPCL_HOST = "pg" ]]; then
#	# borrow env of snpnet
#	libpath="/nfs/data06/ricky/local/lib/snpnet_R4.0.3"
#elif [[ $HPCL_HOST = "obcx" ]]; then
#	libpath="/work/gg57/j29002/.local/lib/R-4.0.3_snpnet2"
#fi

if [ -e ${fout}.lassosum.pipeline.rds ]; then
	echo "rds already exists"
	exit 0
fi

lassosum \
	--out $fout \
	--data $fss --ref.bfile ${gtref} --LDblocks EUR.hg19 \
	--n $nss \
	--chr "CHROM" --pos POS --A1 A1 --A2 A2 --pval P --OR OR \
	--nthreads 1 |&
	tee "${flog}.log"

# ny
script -c "Rscript ./lib/lassosum/src/run_loadrd.R $fname " /dev/null | tee "${flog}_loadrd.log"

date "+%Y-%m-%d %H:%M:%S"

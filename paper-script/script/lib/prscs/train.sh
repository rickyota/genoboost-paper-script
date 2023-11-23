#!/bin/bash
#
# assumed to be called from ./script/nonaddpgs/boosting/pipeline/p1-2_boosting.sh
#

set -eux

program="prscs"

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

if [ $# -ne 5 ]; then
	echo "Wrong argument number."
	exit 1
fi

dout="$1"
ss="$2"
gtref="$3"
nss="$4"
seed="$5"

echo "dout: $dout"
echo "fss: $ss"
echo "gtref: $gtref"
echo "nss: $nss"
echo "seed: $seed"

mkdir -p "$dout"

# TODO
flog="${dout}/${program}.$(date +'%y%m%d-%H%M%S').log"

date "+%Y-%m-%d %H:%M:%S"

pwd

echo "start prscs"

# activate conda before
#python PRScs.py \
script -c "python ./lib/prscs/v210604/PRScs/PRScs.py \
	--ref_dir="./data/prscs/ldblk_1kg_eur" \
	--bim_prefix=$gtref \
	--sst_file=$ss \
	--n_gwas=$nss \
	--out_dir=$dout \
	--seed=$seed" /dev/null |& tee ${flog}

date "+%Y-%m-%d %H:%M:%S"

#!/bin/bash

# first time trial is in ./dataset.12/simu_save/
# at this time, I only used np.random.seed(seed) for the first time,
# and there is no reproducibility.
# Now, I recreate all simu data for reproducibility.

set -eux

dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/conda_activate.sh genetics

#threads=16
#threads=64

ddata="./data/${dataname}/"
ddata_tmp="${ddata}/simu/"
mkdir -p $ddata_tmp

flog="./log/nonaddpgs/dataset/${dataname}/$(basename $0).$(date +'%y%m%d-%H%M%S').log"
mkdir -p "$(dirname $flog)"

#randomn=10

# format
for d in ${ddata_tmp}/*/; do
	[ ! -d $d ] && continue
	dscore="${d}/score-genot/"
	dscore_raw="${dscore}/raw/"

	# TMP
	if [[ $d = *_ncausal-10000* ]]; then

		for f in ${dscore_raw}/*; do
			[ ! -f $f ] && continue
			b=$(basename $f)
			#b=${b/.wgt.score/}
			#mv $f ${dscore}/${b}.score
			mv $f ${dscore}/${b}
		done
	fi
done

# create phenotype
# status col is 'random_0' since plink2 does not allow '-'
python -m projects_py.genetics.dataset.simu \
	--create-phe \
	--dout $ddata_tmp \
	--seed 41 |&
	tee $flog

#source ./script/lib/system/conda_deactivate.sh

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
threads=64

ddata="./data/${dataname}/"
ddata_tmp="${ddata}/simu/"
mkdir -p $ddata_tmp

flog="./log/nonaddpgs/dataset/${dataname}/$(basename $0).$(date +'%y%m%d-%H%M%S').log"
mkdir -p "$(dirname $flog)"

randomn=10

### wgt
python -m projects_py.genetics.dataset.simu \
	--create-wgt \
	--dout $ddata_tmp \
	--fgenot ${ddata}/ukb_imp \
	--genot-format plink2vzs \
	--ffreq ${ddata}/ukb_imp.cv0.afreq \
	--randomn $randomn \
	--seed 31 |&
	tee $flog

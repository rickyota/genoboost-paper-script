#!/bin/bash

set -eux

program="ldpred"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh

#dataname="$1"
#fgenot_ref="$2"

#flog="./log/nonaddpgs/boosting/$(basename $0).$(date +'%y%m%d-%H%M%S').log"
#mkdir -p "$(dirname $flog)"

#dresult="./result/nonaddpgs/${dataname}/${program}/"
#dout="${dresult}/${phe}/train.cv${cvi}/"
#dtmp="${dout}/tmp"
#mkdir -p $dtmp

# unnecessary
#cat ${fgenot_ref}.bim | sed -e "s/^/chrom_/" >${dtmp}1000g.bim
#ln -sf ${fgenot_ref}.fam ${dtmp}/1000g.fam
#ln -sf ${fgenot_ref}.bed ${dtmp}/1000g.bed

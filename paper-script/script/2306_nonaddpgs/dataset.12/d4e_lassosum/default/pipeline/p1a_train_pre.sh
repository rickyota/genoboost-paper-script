#!/bin/bash
# TODO: qsub

set -eux

program="lassosum"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh

dataname="$1"
phes="$2"
cvn="$3"
#fassoc="$2"
#fgenot_ref="$2"

#flog="./log/nonaddpgs/${program}/$(basename $0).$(date +'%y%m%d-%H%M%S').log"
#mkdir -p "$(dirname $flog)"

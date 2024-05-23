#!/bin/bash

# necessary to get nsnv in allstat

set -eux

kind="trans-ancestry_covva"
kind_ori="default"

program="sbayesr"
dataname="dataset.12"

#dscript=$(readlink -f "$(dirname $0)")
# to root
cd "$(dirname $0)/../../../../../"
pwd

source ./script/lib/system/hpc.sh

# shellcheck source=../../../../../data/dataset.12/dataset.config
source "./data/${dataname}/dataset.config"

fpipeline="./script/lib/pipeline/${kind}/p1_train.sh"

phes_wrapper() {
	phes_wrap="$1"
	jobn=$(($(echo "$phes_wrap" | wc -w) * cvn))

	for jobi in $(seq 0 $((jobn - 1))); do
		bash $fpipeline $kind_ori $program $dataname "$phes_wrap" $cvn $jobi
	done

}

phes_wrapper "$phes"

#bash $fpipeline $kind_ori $program $dataname "$phes" $cvn 0

#bash ${dscript}/pipeline/p1_train.sh "$phes" $cvn

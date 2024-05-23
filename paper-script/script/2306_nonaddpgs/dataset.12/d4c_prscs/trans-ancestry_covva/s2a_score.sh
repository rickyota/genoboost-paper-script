#!/bin/bash

set -eux

kind="trans-ancestry_covva"
program="prscs"
dataname="dataset.12"

#dscript=$(readlink -f "$(dirname $0)")
# to root
cd "$(dirname $0)/../../../../../"
pwd

source ./script/lib/system/hpc.sh

# shellcheck source=../../../../../data/dataset.12/dataset.config
source "./data/${dataname}/dataset.config"

fpipeline="./script/lib/pipeline/${kind}/p2a_score.sh"

phes_wrapper() {
	#sex="$1"
	phes_wrap="$1"
	jobn=$(($(echo "$phes_wrap" | wc -w) * cvn))

	for jobi in $(seq 0 $((jobn - 1))); do
		bash $fpipeline $program $dataname "$phes_wrap" $cvn $jobi
	done

}

phes_wrapper "$phes"
#bash $fpipeline $program $dataname "$phes" $cvn 0


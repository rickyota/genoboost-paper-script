#!/bin/bash

set -eux

kind="trans-ancestry_covva"
program="lassosum"
dataname="dataset.12"

#dscript=$(readlink -f "$(dirname $0)")
# to root
cd "$(dirname $0)/../../../../../"
pwd

source ./script/lib/system/hpc.sh

# shellcheck source=../../../../../data/dataset.12/dataset.config
source "./data/${dataname}/dataset.config"

fpipeline="./script/lib/pipeline/${kind}/p2d_score-withcov.sh"

mid_path="${kind}"

phes_wrapper() {
	#sex="$1"
	phes_wrap="$1"
	jobn=$(($(echo "$phes_wrap" | wc -w) * cvn))

	if [ $HPCL_HOST = "pg" ]; then
		#qsub -t 1:$jobn ${dscript}/pipeline/p2d_score-withcov.sh "$phes_wrap" $cvn
		:
	elif [ $HPCL_HOST = "obcx" ]; then
		:
	elif [ $HPCL_HOST = "mdx" ]; then
		sbatch -a 0-$((jobn - 1)) $fpipeline $program $dataname "$mid_path" "$phes_wrap" $cvn
		#sbatch -a 0-$((jobn - 1)) $fpipeline $program $dataname "$phes_wrap" $cvn
		#sbatch -a 0-$((jobn - 1)) ${dscript}/pipeline/p2c_score-withcov_tr.sh $sex "$phes_wrap" $cvn
	fi
}

phes_wrapper "$phes"

#bash $fpipeline $program $dataname "$phes_bothsex" $cvn 0

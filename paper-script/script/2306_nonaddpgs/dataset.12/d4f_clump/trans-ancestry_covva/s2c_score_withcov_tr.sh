#!/bin/bash

set -eux

kind="trans-ancestry_covva"
program="clump"
dataname="dataset.12"

dataname="dataset.12"

#dscript=$(readlink -f "$(dirname $0)")
# to root
cd "$(dirname $0)/../../../../../"
pwd

source ./script/lib/system/hpc.sh

# shellcheck source=../../../../../data/dataset.12/dataset.config
source "./data/${dataname}/dataset.config"

fpipeline="./script/lib/pipeline/${kind}/p2c_score-withcov_tr.sh"

mid_path="${kind}"

phes_wrapper() {
	sex="$1"
	phes_wrap="$2"
	jobn=$(($(echo "$phes_wrap" | wc -w) * cvn))

	if [ $HPCL_HOST = "pg" ]; then
		:
	elif [ $HPCL_HOST = "obcx" ]; then
		:
	elif [ $HPCL_HOST = "mdx" ]; then
		sbatch -a 0-$((jobn - 1)) $fpipeline $program $dataname "$mid_path" $sex "$phes_wrap" $cvn
	fi
}

phes_wrapper both "$phes_bothsex"
phes_wrapper female "$phes_female"

#bash $fpipeline $program $dataname both "$phes_bothsex" $cvn 0

#phes_wrapper() {
#	sex="$1"
#	phes_wrap="$2"
#	jobn=$(($(echo "$phes_wrap" | wc -w) * cvn))
#
#	if [ $HPCL_HOST = "pg" ]; then
#		qsub -t 1:$jobn ${dscript}/pipeline/p2c_score-withcov_tr.sh $sex "$phes_wrap" $cvn
#	elif [ $HPCL_HOST = "obcx" ]; then
#		:
#	elif [ $HPCL_HOST = "mdx" ]; then
#		sbatch -a 0-$((jobn - 1)) ${dscript}/pipeline/p2c_score-withcov_tr.sh $sex "$phes_wrap" $cvn
#	fi
#}
#
#
#phes_wrapper both "$phes_bothsex"
#phes_wrapper female "$phes_female"
#
##bash ${dscript}/pipeline/p2c_score-withcov_tr.sh both "$phes_bothsex" $cvn 0

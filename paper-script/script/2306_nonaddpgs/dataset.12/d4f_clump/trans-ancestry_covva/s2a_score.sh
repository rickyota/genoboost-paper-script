#!/bin/bash

set -eux

dataname="dataset.12"

dscript=$(readlink -f "$(dirname $0)")
# to root
cd "$(dirname $0)/../../../../../"
pwd

source ./script/lib/system/hpc.sh

# shellcheck source=../../../../../data/dataset.12/dataset.config
source "./data/${dataname}/dataset.config"

bash ${dscript}/pipeline/p2a_score.sh "$phes" $cvn

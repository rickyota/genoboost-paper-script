#!/bin/bash
#
# TMP --use-snv-pos
#

set -eux

echo "Start $(readlink -f $0)"
hostname
if [ -v SGE_O_WORKDIR ]; then
    # assume in the root dir
    :
else
    # to root
    cd "$(dirname $0)/../../../"
fi
pwd

source ./script/lib/system/hpc.sh

if [ $# -ne 7 ]; then
    echo "Wrong argument number."
    exit 1
fi

dout="$1"
dtr="$2"
fgenot="$3"
fphe="$4"
lrs="$5"
iters="$6"
threads="$7"

mkdir -p "$dout"

set +ux
if [ $HPCL_HOST = "obcx" ]; then
    module purge
    source ./script/lib/system/module.sh load gcc/7.5.0
fi
set -ux

boosting_wrapper() {
    export RUST_BACKTRACE=full
    ./projects_rust/target/release/boosting_res \
        score \
        --verbose \
        --threads $threads \
        --dir-score $dout \
        --dir-wgt $dtr \
        --file-genot $fgenot \
        --genot-format "plink2-vzs" \
        --file-phe $fphe \
        --learning-rates "$lrs" \
        --iters "$iters" \
        --use-snv-pos \
        "$@"
}

# TODO
# use ukb_imp_pcname.phe
boosting_wrapper --cov sex,age,PC1-PC10
#boosting_wrapper --cov sex,age,PC_1-PC_10

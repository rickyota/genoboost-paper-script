#!/bin/bash
#
# assumed to be called from ./script/nonaddpgs/boosting/pipeline/p1-2_boosting.sh
#

set -eux

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh

if [ $# -ne 11 ]; then
    echo "Wrong argument number."
    exit 1
fi

phe="$1"
cvi="$2"
lrs="$3"
itern="$4"
method="$5"
sex="$6"
dout="$7"
fgenot="$8"
fphe="$9"
fcv="${10}"
threads="${11}"

set +ux
if [ $HPCL_HOST = "obcx" ]; then
    module purge
    source ./script/lib/system/module.sh load gcc/7.5.0
fi
set -ux

mkdir -p "$dout"

# example of $method
# "logitnomissing_effeps-lims2gmodeloverkeepsignprop-100-4-4-0.8_batch-fix-50-10-50_test42"
# -> test42 will be ignroed
# "logitadd_batchfix100-10-100"

# RENEW when option is added
# now adjloss is default
#options="adjloss"
options=""
optiond="effeps batch covway"

IFS='_' read -r -a args_ar <<<"${method}"
args=""
for arg_index in "${!args_ar[@]}"; do
    arg=${args_ar[$arg_index]}
    #echo "$index ${array[index]}"

    isin=false
    if [ $arg_index -eq 0 ]; then
        args="$args --boost-type $arg"
        isin=true
    fi

    # switch option
    for opt in $options; do
        if [ $arg == "${opt}" ]; then
            args="$args --${opt}"
            isin=true
        fi
    done

    # keyword option
    for opt in $optiond; do
        if [[ $arg == ${opt}-* ]]; then
            args="$args --${opt} ${arg/${opt}-/}"
            isin=true
        fi
    done

    if [ ! $isin ]; then
        echo "unused argument: $arg"
    fi

done

echo "args: $args"

phes_wrapper() {
    export RUST_BACKTRACE=full
    ./projects_rust/target/release/boosting_res \
        train \
        --verbose \
        --threads $threads \
        --dir $dout \
        --file-genot $fgenot \
        --genot-format "plink2-vzs" \
        --file-phe $fphe \
        --phe $phe \
        --learning-rates $lrs \
        --iter-snv $itern \
        $args \
        "$@"
}

if [ $sex == "female" ]; then
    phes_wrapper \
        --file-sample <(awk -v cv=${cvi} '($2==0 && $(3+cv) == "tr"){ print $1 }' $fcv) \
        --cov age,PC_1-PC_10
elif [ $sex == "both" ]; then
    phes_wrapper \
        --file-sample <(awk -v cv=${cvi} '($(3+cv) == "tr"){ print $1 }' $fcv) \
        --cov sex,age,PC_1-PC_10
fi

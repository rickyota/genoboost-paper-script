#!/bin/bash

set -eux

dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/conda_activate.sh genetics

dliftover="$1"
fhm3="$2"

ddata="./data/${dataname}/"
#ddata_tmp="${ddata}/tmp_$(basename $(basename $0) .sh)/"
#mkdir -p $ddata_tmp

# done before
# TODO: result should be in ./genetics
fhm3_lift=${fhm3}.liftover
flift_chain="${dliftover}hg18ToHg19.over.chain.gz"
cat ${fhm3}.map | awk ' { print "chr"$1"\t"$4"\t"$4+1"\t"$2 } ' >${fhm3_lift}_hg18.bed
# liftOver only run on z*,g*
${dliftover}liftOver ${fhm3_lift}_hg18.bed $flift_chain ${fhm3_lift}.bed ${fhm3_lift}.unmap.bed

fhm3_lift_bed="${fhm3}.liftover.bed"

snvs_valid_info="${ddata}ukb_imp.valid_info.snvs"
# TOFIX: somehow, rs1377872_G_A is duplicated in snvs_include.
snvs_valid_info_hm3="${ddata}ukb_imp.valid_info_hm3.snvs"

python -m projects_py.genetics.dataset.qc_snv \
    --join_hm3 \
    --fout $snvs_valid_info_hm3 \
    --fhm3 $fhm3_lift_bed \
    --fsnv $snvs_valid_info

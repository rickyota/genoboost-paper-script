#!/bin/bash

set -eux

dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh

#prop_test="$2"
cvn="$1"
#seed="$4"

ddata="./data/${dataname}/"

# all samples are test
# TODO: use cvn
cat ${ddata}/ukb_imp_noneur.psam |
	awk '{ if (NR==1){ print "IID\tsex\tcv0\tcv1\tcv2\tcv3\tcv4" }else{ print $1"\t"2-$2"\tts\tts\tts\tts\tts" }}' \
		>${ddata}/ukb_imp_noneur.cv

#ddata_tmp="${ddata}/tmp/"
#mkdir -p $ddata_tmp

#python -m projects_py.genetics.dataset.cross_vali \
#	--cv \
#	--fout ${ddata}/ukb_imp.cv \
#	--fgenot ${ddata}/ukb_imp \
#	--genot-format plink2 \
#	--prop_test $prop_test \
#	--cvn $cvn \
#	--seed $seed
#
#source ./script/lib/system/conda_deactivate.sh

#!/bin/bash

set -eux

dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/conda_activate.sh genetics

#prop_test="$2"
cvn="$1"
seed="$2"

ddata="./data/${dataname}/"

# all samples are test
# TODO: use cvn
#cat ${ddata}/ukb_imp_noneur.psam |
#	awk '{ if (NR==1){ print "IID\tsex\tcv0\tcv1\tcv2\tcv3\tcv4" }else{ print $1"\t"2-$2"\tts\tts\tts\tts\tts" }}' \
#		>${ddata}/ukb_imp_noneur.cv

ddata_tmp="${ddata}/tmp/"
mkdir -p $ddata_tmp

# only set one kind of va to make all test samples common

for anc in afr sas eas; do
	## FIXME: python: make cvn=1 work
	## TMP: output va and ts only
	python -m projects_py.genetics.dataset.cross_vali \
		--cv \
		--fout ${ddata_tmp}/ukb_imp_trans-ancestry_vats_cv0.${anc}.cv \
		--fgenot ${ddata}/ukb_imp_trans-ancestry \
		--genot-format plink2 \
		--fsamples ${ddata}/qc.trans-ancestry.${anc}.sample \
		--prop-test 0.8 \
		--cvn 2 \
		--seed $seed

	# TODO: use cvn
	# only use cv0 and copy to cv1-4
	# make tr->va
	cat ${ddata_tmp}/ukb_imp_trans-ancestry_vats_cv0.${anc}.cv |
		awk '{ if (NR==1){ print "IID\tsex\tcv0\tcv1\tcv2\tcv3\tcv4" }else{ print $1"\t"$2"\t"$3"\t"$3"\t"$3"\t"$3"\t"$3 }}' |
		sed -e "s/tr/va/g" \
			>${ddata}/ukb_imp_trans-ancestry_vats.${anc}.cv
done

# unnecessary?
awk 'FNR>1 || NR==1' ${ddata}/ukb_imp_trans-ancestry_vats.{afr,sas,eas}.cv >${ddata}/ukb_imp_trans-ancestry_vats.cv

## TODO: use cvn
## only use cv0 and copy to cv1-4
## make tr->va
#cat ${ddata_tmp}/ukb_imp_trans-ancestry_vats.cv |
#	awk '{ if (NR==1){ print "IID\tsex\tcv0\tcv1\tcv2\tcv3\tcv4" }else{ print $1"\t"$2"\t"$3"\t"$3"\t"$3"\t"$3"\t"$3 }}' |
#	sed -e "s/tr/va/g" \
#		>${ddata}/ukb_imp_trans-ancestry_vats.cv

#tr 'tr' 'va' \

#!/bin/bash

set -eux

echo "Start $(readlink -f $0)"
hostname
if [ -v SGE_O_WORKDIR ]; then
	# assume you are in the root
	:
else
	# to root
	cd "$(dirname $0)/../../../"
fi
pwd

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh

if [ $# -ne 11 ]; then
	echo "Wrong argument number."
	exit 1
fi

batch_id="$1"
gmodel="$2"
cvi="$3"
dout="$4"
fgenot="$5"
fphe="$6"
fcv="$7"
sex="$8"
fbatch="$9"
phes="${10}"
threads="${11}"

declare -A gmodel_plink=(
	["add"]=""
	["rec"]="recessive"
	["dom"]="dominant"
	["hetonly"]="hetonly"
)

gm_p=${gmodel_plink[$gmodel]}

source ./script/lib/system/module.sh load plink2

if [ $gmodel == "add" ]; then
	fout="${dout}ukb_imp.cv${cvi}.batch${batch_id}"
else
	fout="${dout}ukb_imp.cv${cvi}.batch${batch_id}.${gmodel}"
fi

plink_glm_wrapper() {
	# failed: did not converge
	#--glm zs no-x-sex log10 hide-covar skip-invalid-pheno cc-residualize single-prec-cc firth-fallback \
	#--glm $gm_p zs no-x-sex log10 hide-covar skip-invalid-pheno firth-fallback \
	plink2 \
		--glm $gm_p zs omit-ref no-x-sex hide-covar skip-invalid-pheno firth-fallback \
		--ci 0.95 \
		--out $fout \
		--pfile $fgenot vzs \
		--chr 1-22,X,XY,Y,MT \
		--extract <(awk -v batch_id=${batch_id} '($NF == batch_id){ print $3 }' $fbatch) \
		--covar ${fphe} \
		--covar-variance-standardize \
		--pheno-quantile-normalize \
		--1 \
		--threads $threads \
		--memory 16000 \
		--pheno ${fphe} \
		--pheno-name ${phes} \
		"$@"
}

if [ $sex == "both" ]; then
	plink_glm_wrapper \
		--keep <(awk -v cv=${cvi} '($(3+cv) == "tr"){ print $1 }' $fcv) \
		--covar-name sex,age,PC_1-PC_10
elif [ $sex == "female" ]; then
	plink_glm_wrapper \
		--keep <(awk -v cv=${cvi} '($2==0 && $(3+cv) == "tr"){ print $1 }' $fcv) \
		--covar-name age,PC_1-PC_10
else
	echo "Unknown sex: $sex"
fi

source ./script/lib/system/module.sh unload plink2

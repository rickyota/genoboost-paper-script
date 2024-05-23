#!/bin/bash

# boosting and snpnet is not using covva regression

set -eux

id="2_snpboost"
kind="trans-ancestry_covva"
program="paper"
dataname="dataset.12"
#model="noneur"

#dscript=$(readlink -f "$(dirname $0)")
# to root
cd "$(dirname $0)/../../../../../"

# shellcheck source=../../../../../data/dataset.12/dataset.config
source "./data/${dataname}/dataset.config"

runs="acc-trans-ancestry-covva-boost"
#runs="acc-noneur-covva-boost"

#methods="boosting_integrate boosting_nonadd boosting_add snpnet lassosum ldpred prscs sbayesr clump \
#boosting_integrate-noneur boosting_nonadd-noneur boosting_add-noneur \
#snpnet_noneur lassosum_noneur ldpred_noneur prscs_noneur sbayesr_noneur clump_noneur"

model_boost_nonadd="logitnomissing_covway-first_effeps-lims2gmodeloverkeepsignprop-100-4-4-0.8_batch-fix-50-10-50_gwasmethods"
model_boost_add="logitadd_covway-first_batch-fix-50-10-50_gwasmethods"
#model_boost_nonadd_ex_chr6="logitnomissing_covway-first_effeps-lims2gmodeloverkeepsignprop-100-4-4-0.8_batch-fix-50-10-50_ex-chr6_230507"
#model_boost_add_ex_chr6="logitadd_covway-first_batch-fix-50-10-50_ex-chr6_230507"
model_boost_loss="logitnomissing_covway-first_effeps-lims2gmodeloverkeepsignprop-100-4-4-0.8_batch10000000_loss_gwasmethods"

# TMP
#phes="bc"

# for order
phes="ra psr gout ibd atm acd ad af bc cc cad t2d"
#phes="ra psr"
#phes="ra atm psr gout ibd bc acd ad cc af cad t2d her ingher ami mch mi hae mbe"
#phes="t2d cad bc af ibd ra gout ad acd cc atm psr her ingher ami mch mi hae mbe"

#bash ${dscript}/pipeline/p1_paper.sh $dataname $cvn $phes $phes_female $methods \
#	$model_boost_nonadd $model_boost_add \
#	$model_boost_loss $fgene_exon

threads=16

para_cand="cand13"
para_com="com1"
para_regon="regon4"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/conda_activate.sh genetics
source ./script/lib/system/python_threads.sh $threads

ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/"

fcv="${ddata}/ukb_imp.cv"
fphe="${ddata}/ukb_imp.phe"

fgene_exon="./data/ensembl/Homo_sapiens.GRCh37.87.gtf.ensembl.gene_exon.protein_coding.tsv"

dout="${dresult}/paper/${kind}/${model_boost_nonadd}_${model_boost_add}_${para_cand}_${para_regon}_${id}/"
mkdir -p $dout

flog="./log/nonaddpgs/${dataname}/${program}/${kind}/${model_boost_nonadd}_${model_boost_add}_${para_cand}_${para_regon}_${id}/paper.$(date +'%y%m%d-%H%M%S').log"
mkdir -p "$(dirname $flog)"

#for anc in eas; do
#for anc in afr sas eas; do
for anc in sas; do
	methods="boosting_integrate boosting_nonadd boosting_add snpboost snpnet lassosum ldpred prscs sbayesr clump \
	boosting_integrate-${anc} boosting_nonadd-${anc} boosting_add-${anc} \
	snpboost_${anc} snpnet_${anc} lassosum_${anc} ldpred_${anc} prscs_${anc} sbayesr_${anc} clump_${anc}"
	methods_type="boosting_integrate-${anc} boosting_nonadd-${anc} boosting_add-${anc} \
	snpboost_${anc} snpnet_${anc} lassosum_${anc} ldpred_${anc} prscs_${anc} sbayesr_${anc} clump_${anc}"

	python -m projects_py.boosting.main.paper \
		--format $program \
		--runs $runs \
		--dout $dout \
		--dresult $dresult \
		--ddata $ddata \
		--method-mid-path "boosting_nonadd" default $model_boost_nonadd \
		--method-mid-path "boosting_add" default $model_boost_add \
		--method-mid-path "boosting_nonadd-${anc}" trans-ancestry_covva $model_boost_nonadd ${anc} \
		--method-mid-path "boosting_add-${anc}" trans-ancestry_covva $model_boost_add ${anc} \
		--method-mid-path "boosting_loss" default $model_boost_loss \
		--method-mid-path "snpboost" monitor \
		--method-mid-path "snpboost_${anc}" trans-ancestry_covva ${anc} \
		--method-mid-path "snpnet_${anc}" trans-ancestry_covva ${anc} \
		--method-mid-path "lassosum_${anc}" trans-ancestry_covva ${anc} \
		--method-mid-path "ldpred_${anc}" trans-ancestry_covva ${anc} \
		--method-mid-path "prscs_${anc}" trans-ancestry_covva ${anc} \
		--method-mid-path "sbayesr_${anc}" trans-ancestry_covva ${anc} \
		--method-mid-path "clump_${anc}" trans-ancestry_covva ${anc} \
		--model-integrate "boosting_integrate-${anc}" "boosting_add-${anc} boosting_nonadd-${anc}" \
		--methods $methods \
		--methods-type ${anc} $methods_type \
		--fcv $fcv \
		--cvn $cvn \
		--fpheno $fphe \
		--phes $phes --phes-female $phes_female \
		--para-cand $para_cand --para-com $para_com --para-regon $para_regon \
		--fgene-exon $fgene_exon |&
		tee $flog
done

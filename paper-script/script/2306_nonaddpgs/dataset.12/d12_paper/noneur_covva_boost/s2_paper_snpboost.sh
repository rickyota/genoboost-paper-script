#!/bin/bash

# boosting and snpnet is not using covva regression

set -eux

id="2_snpboost"
kind="noneur_covva_boost"
program="paper"
dataname="dataset.12"

# to root
cd "$(dirname $0)/../../../../../"

# shellcheck source=../../../../../data/dataset.12/dataset.config
source "./data/${dataname}/dataset.config"

runs="acc-noneur-covva-boost"

methods="boosting_integrate boosting_nonadd boosting_add snpboost snpnet lassosum ldpred prscs sbayesr clump \
boosting_integrate-noneur boosting_nonadd-noneur boosting_add-noneur \
snpboost_noneur snpnet_noneur lassosum_noneur ldpred_noneur prscs_noneur sbayesr_noneur clump_noneur"
methods_type="boosting_integrate-noneur boosting_nonadd-noneur boosting_add-noneur \
snpnet_noneur lassosum_noneur ldpred_noneur prscs_noneur sbayesr_noneur clump_noneur"

model_boost_nonadd="logitnomissing_covway-first_effeps-lims2gmodeloverkeepsignprop-100-4-4-0.8_batch-fix-50-10-50_gwasmethods"
model_boost_add="logitadd_covway-first_batch-fix-50-10-50_gwasmethods"
model_boost_loss="logitnomissing_covway-first_effeps-lims2gmodeloverkeepsignprop-100-4-4-0.8_batch10000000_loss_gwasmethods"

# for order
phes="ra psr gout ibd atm acd ad af bc cc cad t2d"

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

python -m projects_py.boosting.main.paper \
	--runs $runs \
	--dout $dout \
	--dresult $dresult \
	--ddata $ddata \
	--method-mid-path "boosting_nonadd" default $model_boost_nonadd \
	--method-mid-path "boosting_add" default $model_boost_add \
	--method-mid-path "boosting_nonadd-noneur" noneur_covva $model_boost_nonadd \
	--method-mid-path "boosting_add-noneur" noneur_covva $model_boost_add \
	--method-mid-path "boosting_loss" default $model_boost_loss \
	--method-mid-path "snpboost" monitor \
	--method-mid-path "snpboost_noneur" noneur_covva \
	--method-mid-path "snpnet_noneur" noneur_covva \
	--method-mid-path "lassosum_noneur" noneur_covva \
	--method-mid-path "ldpred_noneur" noneur_covva \
	--method-mid-path "prscs_noneur" noneur_covva \
	--method-mid-path "sbayesr_noneur" noneur_covva \
	--method-mid-path "clump_noneur" noneur_covva \
	--model-integrate "boosting_integrate-noneur" "boosting_add-noneur boosting_nonadd-noneur" \
	--methods $methods \
	--methods-type noneur $methods_type \
	--fcv $fcv \
	--cvn $cvn \
	--fpheno $fphe \
	--phes $phes --phes-female $phes_female \
	--para-cand $para_cand --para-com $para_com --para-regon $para_regon \
	--fgene-exon $fgene_exon |&
	tee $flog

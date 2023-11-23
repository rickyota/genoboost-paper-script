#!/bin/bash

# all relative path should be from root (./genetics/)

set -eux

id="1"
kind="ex-chr6"
program="paper"
dataname="dataset.12"

# to root
cd "$(dirname $0)/../../../../../"
pwd

# shellcheck source=../../../../../data/dataset.12/dataset.config
source "./data/${dataname}/dataset.config"

runs="acc"

methods="boosting_integrate boosting_nonadd boosting_add ldpred \
boosting_integrate-ex-chr6 boosting_nonadd-ex-chr6 boosting_add-ex-chr6 ldpred_ex-chr6"
methods_type="boosting_integrate-ex-chr6 boosting_nonadd-ex-chr6 boosting_add-ex-chr6 ldpred_ex-chr6"
model_boost_nonadd="logitnomissing_covway-first_effeps-lims2gmodeloverkeepsignprop-100-4-4-0.8_batch-fix-50-10-50_gwasmethods"
model_boost_add="logitadd_covway-first_batch-fix-50-10-50_gwasmethods"
model_boost_nonadd_ex_chr6="logitnomissing_covway-first_effeps-lims2gmodeloverkeepsignprop-100-4-4-0.8_batch-fix-50-10-50_ex-chr6_230507"
model_boost_add_ex_chr6="logitadd_covway-first_batch-fix-50-10-50_ex-chr6_230507"
model_boost_loss="logitnomissing_covway-first_effeps-lims2gmodeloverkeepsignprop-100-4-4-0.8_batch10000000_loss_gwasmethods"

# for order
phes="ra psr"

fgene_exon="./data/ensembl/Homo_sapiens.GRCh37.87.gtf.ensembl.gene_exon.protein_coding.tsv"

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

dout="${dresult}/${program}/${kind}/${model_boost_nonadd}_${model_boost_add}_${para_cand}_${para_regon}_${id}/"
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
	--method-mid-path "boosting_nonadd-ex-chr6" ex-chr6 $model_boost_nonadd_ex_chr6 \
	--method-mid-path "boosting_add-ex-chr6" ex-chr6 $model_boost_add_ex_chr6 \
	--method-mid-path "ldpred_ex-chr6" ex-chr6-gwasmethods \
	--method-mid-path "boosting_loss" default $model_boost_loss \
	--model-integrate "boosting_integrate-ex-chr6" "boosting_add-ex-chr6 boosting_nonadd-ex-chr6" \
	--methods-type ex-chr6 $methods_type \
	--methods $methods \
	--fcv $fcv \
	--cvn $cvn \
	--fpheno $fphe \
	--phes $phes --phes-female $phes_female \
	--para-cand $para_cand --para-com $para_com --para-regon $para_regon \
	--fgene-exon $fgene_exon |& tee $flog

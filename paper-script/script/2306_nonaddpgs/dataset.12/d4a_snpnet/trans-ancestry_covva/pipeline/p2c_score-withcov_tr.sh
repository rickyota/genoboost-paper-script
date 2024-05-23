#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -q all.q
#$ -l hostname=a02|a03|z01|z02|g01|g02

# slurm
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=30MB
#SBATCH --partition=all

set -eux

kind="trans-ancestry_covva"
#mode="noneur"
program="snpnet"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh
source ./script/lib/system/module_init.sh
source ./script/lib/system/conda_activate.sh genetics

if [ -v SGE_TASK_ID ]; then
	job_id=$((SGE_TASK_ID - 1))
	threads=$NSLOTS
elif [ -v SLURM_ARRAY_TASK_ID ]; then
	job_id=$SLURM_ARRAY_TASK_ID
	threads=$SLURM_CPUS_PER_TASK
else
	job_id=$4
	threads=16
fi

source ./script/lib/system/python_threads.sh $threads

sex="$1" # male, female, both
phes="$2"
cvn="$3"

cvi=$((job_id % cvn))
phei=$((job_id / cvn))
phes_ar=(${phes})
phe=${phes_ar[$phei]}

ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/"

#fcv="${ddata}ukb_imp_noneur_vats.cv"
fphe="${ddata}ukb_imp_trans-ancestry_pcname.phe"
#fphe="${ddata}ukb_imp_noneur_pcname.phe"

if [ $sex == "female" ]; then
	cov_names="age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
elif [ $sex == "both" ]; then
	cov_names="sex,age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
fi

flog="./log/nonaddpgs/${dataname}/${program}/${kind}/${phe}/score.cv${cvi}/score-withcov.$(date +'%y%m%d-%H%M%S').log"
mkdir -p "$(dirname $flog)"

for anc in afr sas eas; do
	fcv="${ddata}ukb_imp_trans-ancestry_vats.${anc}.cv"

	python -m projects_py.boosting.main.withcov \
		--train \
		--regon va \
		--method $program \
		--mid-path $kind $anc \
		--dresult $dresult \
		--fcv $fcv --cvi $cvi \
		--fpheno $fphe --phe $phe \
		--cov-names $cov_names \
		--sex $sex |& tee $flog
done

#bash ./script/${dataname}/lib/${program}/score-withcov.sh \
# $phe $cvi $dresult $fphe $fcv $kind  |& tee $flog

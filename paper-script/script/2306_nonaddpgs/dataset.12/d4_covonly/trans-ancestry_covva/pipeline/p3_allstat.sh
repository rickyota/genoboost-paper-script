#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -q all.q
#$ -l hostname=z02
##$ -l hostname=a02|a03|z01|z02|g01|g02

# slurm
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=30MB
#SBATCH --partition=all

set -eux

kind="trans-ancestry_covva"
program="covonly"
dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

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

sex="$1"
phes="$2"
cvn="$3"

regon='va'
regon_covonly='va'

phei=$job_id
phes_ar=(${phes})
phe=${phes_ar[$phei]}

ddata="./data/${dataname}/"
dresult="./result/nonaddpgs/${dataname}/"
mkdir -p $dresult

fphe="${ddata}ukb_imp_trans-ancestry.phe"
#fcv="${ddata}ukb_imp_noneur_vats.cv"
#fphe="${ddata}ukb_imp_noneur.phe"

flog="./log/nonaddpgs/${dataname}/${program}/${kind}/${phe}/allstat.$(date +'%y%m%d-%H%M%S').log"
mkdir -p "$(dirname $flog)"

for anc in afr sas eas; do
    fcv="${ddata}ukb_imp_trans-ancestry_vats.${anc}.cv"

    python -m projects_py.boosting.main.allstat \
        --dresult $dresult \
        --method $program \
        --mid-path $kind $anc \
        --mid-path-covonly trans-ancestry_covva $anc \
        --regon $regon \
        --regon_covonly $regon_covonly \
        --fcv $fcv --cvn $cvn \
        --fpheno $fphe --phe $phe \
        --sex $sex |& tee $flog
done

## TODO: afr_sas_eas anc
# create anc="afr_sas_eas" dir and put merged samples there

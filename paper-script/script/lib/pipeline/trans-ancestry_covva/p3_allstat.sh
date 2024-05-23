#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 2
#$ -q all.q
#$ -l hostname=a02|a03|z01|z02|g01|g02

# slurm
#SBATCH --cpus-per-task=32
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=30MB
#SBATCH --partition=all

set -eux

#kind="tranancestry_covva"
#program="boosting"
#dataname="dataset.12"

if [ $# -ne 8 ] && [ $# -ne 9 ]; then
    echo "Wrong argument number."
    exit 1
fi

program=$1
dataname=$2
mid_path=$3 # ${kind}/${model}

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
    job_id=$9
    threads=32
fi

source ./script/lib/system/python_threads.sh $threads

sex="$4"
phes="$5"
cvn="$6"
#model="$4"

regon="$7"
regon_covonly="$8"
#regon='va'
#regon_covonly='va'

phei=$job_id
phes_ar=(${phes})
phe=${phes_ar[$phei]}

ddata="./data/${dataname}/"
#dresult="./result/nonaddpgs/${dataname}/${program}/"
dresult="./result/nonaddpgs/${dataname}/"

#fcv="${ddata}ukb_imp_noneur_vats.cv"
#fcv="${ddata}/ukb_imp_noneur.cv"
# depends on method?
fphe="${ddata}ukb_imp_trans-ancestry_pcname.phe"
#fphe="${ddata}ukb_imp_noneur_pcname.phe"
#fphe="${ddata}/ukb_imp_noneur.phe"

#fcv="${ddata}/ukb_imp.cv"
#fphe="${ddata}/ukb_imp.phe"

flog="./log/nonaddpgs/${dataname}/${program}/${mid_path}/${phe}/allstat.$(date +'%y%m%d-%H%M%S').log"
mkdir -p "$(dirname $flog)"

#for anc in eas; do
for anc in afr sas eas; do
    fcv="${ddata}ukb_imp_trans-ancestry_vats.${anc}.cv"

    python -m projects_py.boosting.main.allstat \
        --dresult $dresult \
        --method $program \
        --mid-path $mid_path $anc \
        --mid-path-covonly trans-ancestry_covva $anc \
        --regon $regon \
        --regon_covonly $regon_covonly \
        --fcv $fcv --cvn $cvn \
        --fpheno $fphe --phe $phe \
        --sex $sex |& tee $flog

    #--mid-path $kind $model $anc \

done

## TODO: afr_sas_eas anc
# create anc="afr_sas_eas" dir and put merged samples there

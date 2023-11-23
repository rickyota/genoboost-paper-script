#!/bin/bash

# copy genot file to SSD in pg

set -eux

dataname="dataset.12"

hostname
pwd
echo "Start $(readlink -f $0)"

source ./script/lib/system/hpc.sh



if [ $HPCL_HOST != "pg" ]; then
    echo "run only on pg"
    exit 0
fi

#nodes="z01 z02 a02 a03 g01"
#nodes="g01"
nodes="g02"

ddata="./data/${dataname}/"

# create dir
for node in $nodes; do
    echo $node
    ssh $node <<EOF
        mkdir -p /grid2/ricky/data/${dataname}/
EOF
done

cd /nfs/data06/ricky/code/genetics

for node in $nodes; do
    echo $node
    ddata="${node}:/grid2/ricky/data/${dataname}/"
    rsync -ahvz ./data/${dataname}/ukb_imp_noneur.{pgen,psam,pvar.zst} $ddata
    rsync -ahvz ./data/${dataname}/ukb_imp.{pgen,psam,pvar.zst} $ddata
    #rsync -ahvz ./data/${dataname}/ukb_imp.{bed,bim,fam} $ddata
    #  --copy-links
done

# do not sync .phe etc. since it could change a lot later

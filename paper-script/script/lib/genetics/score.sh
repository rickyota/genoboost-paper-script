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

if [ $# -ne 4 ]; then
	echo "Wrong argument number."
	exit 1
fi

dout="$1"
dwgt="$2"
fgenot="$3"
threads="$4"

set +ux
if [ $HPCL_HOST = "obcx" ]; then
	module purge
	source ./script/lib/system/module.sh load gcc/7.5.0
fi
set -ux

export RUST_BACKTRACE=full

./projects_rust/target/release/genetics_res \
	score \
	--verbose \
	--threads $threads \
	--dir-score ${dout} \
	--file-genot ${fgenot} \
	--genot-format "plink2-vzs" \
	--dir-wgt ${dwgt}

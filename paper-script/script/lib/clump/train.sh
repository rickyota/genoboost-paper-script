#!/bin/bash
#
# assumed to be called from ./script/nonaddpgs/boosting/pipeline/p1-2_boosting.sh
#

set -eux

program="clump"

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
module load plink

if [ $# -ne 11 ]; then
	echo "Wrong argument number."
	exit 1
fi

phe="$1"
cvi="$2"
dout="$3"
fgenot="$4"
fphe="$5"
fassoc="$6"
fcv="$7"
sex="$8"
threads="$9"
# usually run with max of ps and extract later
clump_ps="${10}"
clump_r2s="${11}"

echo "dout: $dout"

mkdir -p "$dout"

clump_kb=250

# TODO
flog="${dout}/${program}.$(date +'%y%m%d-%H%M%S').log"

for clump_p in $clump_ps; do
	# do not use p2 -> any small number
	clump_p2=$clump_p

	for clump_r2 in $clump_r2s; do
		echo "clump_p,clump_r2: $clump_p , $clump_r2"

		#fi

		fout="${dout}/clump.p${clump_p}.r${clump_r2}"

		if [ -e ${fout}.clumped ]; then
			echo ".clumped exists"
			continue
		fi

		train_wrapper() {
			plink --out $fout \
				--threads ${threads} \
				--bfile ${fgenot} \
				--clump $fassoc \
				--clump-snp-field ID --clump-p1 $clump_p --clump-p2 $clump_p2 --clump-r2 $clump_r2 --clump-kb $clump_kb
		}

		if [ $sex == "both" ]; then
			train_wrapper \
				--keep <(awk -v cv=${cvi} '($(3+cv) == "tr"){ print $1 }' $fcv)
		elif [ $sex == "female" ]; then
			train_wrapper \
				--keep <(awk -v cv=${cvi} '($2==0 && $(3+cv) == "tr"){ print $1 }' $fcv)
		else
			echo "Unknown sex: $sex"
		fi

	done
done

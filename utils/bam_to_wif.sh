#!/usr/bin/env bash

# save a variant file with prefix (1) using a vcf file (2)
function get_variants {
    python get.variants.py $1 $2
}

# convert a bam alignment file (1), using a variant file (2), to a wif file ( for a specific chromosome (3)).
function bam_to_wif {
    python bam_to_wif.py -b $1 -vf $2 -chr $3
}


###PARSER###

if [ $# -eq 0 ]; then
echo "Usage: $0 -p [prefix] -v [vcf file] -b [bam file] -c [chromosome]";
echo "";
echo "[prefix] = The prefix for the variant file";
echo "[vcf file] = Vcf file";
echo "[bam file] = bam alignment file ( need a bai index )";
echo "[chromosome] = The chromosome to analize, default chr1 ";
echo "";
exit 1;
fi

while [[ $# -gt 1 ]]
do
key="$1"

case ${key} in
    -p|--prefix)
    PREFIX="$2"
    shift # past argument
    ;;
    -v|--vcf)
    VCFFILE="$2"
    shift # past argument
    ;;
    -b|--bam)
    BAMFILE="$2"
    shift # past argument
    ;;
    -c|--chromosome)
    CHROMOSOME="$2"
    shift
    ;;
    *)
esac
shift # past argument or value
done
echo prefix = "${PREFIX}"
echo vcf_file = "${VCFFILE}"
echo bam_file = "${BAMFILE}"
echo chromosome = "${CHROMOSOME}"

## MAIN ##

get_variants ${PREFIX} ${VCFFILE}
bam_to_wif ${BAMFILE} ${PREFIX}_${CHROMOSOME}.var ${CHROMOSOME}

#!/usr/bin/env bash

# To execute this script is required a subdirectory structure as described below.:

    # - main_folder
    #   - son
    #       - Tech 1
    #           - alignment.txt ( contains a line with the url of alignment.bam file [required] )
    #           - alignment_index.txt ( contains a line with the url of alignment.bam.bai file [optional] )
    #           - variants.txt ( contains a line with the url of variants.vcf file [required] )
    #       - Tech 2
    #       - Tech 3
    #   - mother
    #   - father
    #       - Tech 1

# script_folder contains the script that will be executed in pipeline

# build variant file using individual (1), chromosome (2), script folder (3) and vcf file (4)
function build_variant {

    if [ ! -f "${1}"_chr"${2}".var ] && [ ! -f "${1}"_"${2}".var ]; then
      {
        if [ ! -z "${4}" ]; then
          var="${4}"
        else
          var="./variants.vcf"
        fi
        python2 "${3}"get.variants.py "${1}" "${var}"
      } || echo 1
    fi

    # Some variants are saved as 'individual'_chr_*, while others such as 'individual'_*

    individual_chr="${1}"_chr"${2}".var

    if [ -f "${1}"_"${2}".var ] ; then

      individual_chr="${1}"_"${2}".var
      echo "${individual_chr}";

    else
        echo 1
    fi
}

# Using individual (1), chromosome (2) and script folder (3) build a sfi file
function build_sfi {

    if [ ! -f "${1}".sfi ]; then
      #
      # get.sfi_pysam uses pysam to fetch content of bam alignment file. If something goes wrong, fallback
      # on "old" script get.sfi.py in pipeline with samtools view
      {
        {
          python "${3}"get.sfi_pysam.py "${2}" alignment.bam > "${1}".sfi
        } ||  # catch
        {
          echo "Fallback. Trying to use samtools and legacy get.sfi.py"
          echo
          samtools view alignment.bam | python2 "${3}"get.sfi.py "${2}" > "${1}".sfi
        } || echo 1
      }
    fi
}

#Using individual (1), chromosome (2) and script folder (3) build matrix with snps info and its transpose
function build_matrix_transpose {

    sanity_check=0
    if [ ! -f "${1}.matrix" ] ; then
      {
        python2 "${3}"get.matrix.py "${2}" "${1}".sfi > "${1}".matrix
      } || sanity_check=1
    fi

    if [ "${sanity_check}" ] && [ ! -f "${1}.transpos" ] ; then
      {
        python2 "${3}"get.transpos.py "${2}" "${1}".matrix > "${1}".transpos
      } || echo 1
    fi

}

# Using individual (1) and script folder build variant statistics
function build_var_stats {

    if [ ! -f "${1}.stats" ] ; then
      {
        python2 "${2}"get.var_stats.py "${1}".transpos > "${1}".stats
      } || echo 1
    fi
}

# Individual (1) data will be deleted. If skip (2) is true, the deletion is automatic.
#
function delete_data {

   echo "Do you want to clear bam, sfi and other chr data? [y/N]"
   read -r input

    if [[ ${input} =~ ^[Yy]$ ]]; then {
      rm -rf alignment.bam alignment.bam.bai variants.vcf "${1}_*" "${1}.sfi" "${1}.matrix" "${1}.transpos"
    }
    fi
}

###PARSER###

if [ $# -eq 0 ]; then
echo "Usage: $0 -i [individual] -m [main folder] -s [script folder] -c [chromosome] -v [vcf file]";
echo "";
echo "[individual] = son/mother/father which data will be analyzed";
echo "[main folder] = the root folder, where son, mother, father subfolders are created";
echo "[script folder] = where the python scripts are located";
echo "[chromosome] = The chromosome to analize, default chr1";
echo "[vcf file] = Specific vcf file with snps position. By default will be downloaded the file with url in variant.txt";
echo "";
exit 1;
fi

while [[ $# -gt 1 ]]
do
key="$1"

case ${key} in
    -i|--individual)
    INDIVIDUAL="$2"
    shift # past argument
    ;;
    -m|--main)
    MAINFOLDER="$2"
    shift # past argument
    ;;
    -s|--script)
    SCRIPTFOLDER="$2"
    shift # past argument
    ;;
    -c|--chromosome)
    CHROMOSOME="$2"
    shift
    ;;
    -v|--vcf-file)
    VCFFILE="$2"
    shift
    ;;
    *)
esac
shift # past argument or value
done
echo individual = "${INDIVIDUAL}"
echo main_folder = "${MAINFOLDER}"
echo script_folder = "${SCRIPTFOLDER}"
echo chromosome = "${CHROMOSOME}"
echo vcf_file = "${VCFFILE}"
echo

## MAIN ##

echo "Pipeline started"
echo

individual_folder="${MAINFOLDER}${INDIVIDUAL}/"

cd "${individual_folder}" || exit

#
# Iterate over individual subfolders (technologies)
#
for path in ${individual_folder}*; do
  [ -d "${path}" ] || continue # if not a directory, skip
  dirname="$(basename "${path}")"

  echo "Entering folder ${dirname}"
  echo

  echo "Do you want to continue with this technology? [y/N]"
  read -r input

  if [[ ${input} =~ ^[Yy]$ ]]; then

    cd "${path}" || exit
    sanity_check=0

    ## Download data ##
    echo "Downloading data"
    echo

    if [ ! -f "alignment.bam" ]; then
      filename='alignment.txt'
      { while read -r url; do
            wget "${url}" -O "alignment.bam" ;
        done < ${filename}
      } || sanity_check=1
    fi

    if [ ! -f "alignment.bam.bai" ] && [ -f "alignment_index.txt" ] ; then
      {
        filename='alignment_index.txt'
        while read -r url; do
          wget "${url}" -O "alignment.bam.bai"
        done < ${filename}
      } || echo "No alignment index file found"

    fi

    if [ -z "${VCFFILE// }" ] && [ ! -f "variants.vcf" ]; then
      filename='variants.txt'
      { while read -r url; do
            wget "${url}" -O "variants.vcf" ;
        done < ${filename}
      } || sanity_check=1
    fi

    if [ ${sanity_check} = 0 ]; then
        echo "Step 1 - Extracting variants from vcf"
        echo

        chr="$(build_variant "${INDIVIDUAL}" "${CHROMOSOME}" "${SCRIPTFOLDER}" "${VCFFILE}")"
        if ! [ "${chr}" = 1 ]; then
            echo "Variants built"
            echo
            echo "Step 2 - Building sfi file"
            echo
            res="$(build_sfi "${INDIVIDUAL}" "${chr}" "${SCRIPTFOLDER}")"
            if [ ! "${res}" = 1 ]; then
                echo "Sfi file built"
                echo
                echo "Step 3 - Building matrix and its transpose"
                echo
                res="$(build_matrix_transpose "${INDIVIDUAL}" "${chr}" "${SCRIPTFOLDER}")"
                if [ ! "${res}" = 1 ]; then
                    echo "Matrix and transpose built"
                    echo
                    echo "Step 4 - Building variants statistics file"
                    echo
                    res="$(build_var_stats "${INDIVIDUAL}" "${SCRIPTFOLDER}")"
                    if [ ! "${res}" = 1 ]; then
                        echo "Variants statistics file created. Now can be deleted all useless files."
                        echo
                        delete_data "${INDIVIDUAL}"
                        echo "Computation finished for ${dirname}"
                        echo
                    else
                        echo "Error during building variants statistics"
                    fi
                else
                    echo "Error during building matrix or transpose"
                fi
            else
                echo "Error during building sfi file"
            fi
        else
            echo "Error during building variant"
        fi
    else
        echo "Error during downloading data"
    fi
  fi
done

echo "Pipeline finished"
echo
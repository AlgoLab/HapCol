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

script_folder="/home/prj_rnabwt/haplotyping/giab/scripts/"
main_folder="/data/AshkenazimTrio/"

echo "Insert individual type: mother - son - father"
read -r individual

individual_folder="${main_folder}${individual}/"

cd "${individual_folder}" || exit

#
# Iterate over individual subfolders (technologies)
#
for path in ${individual_folder}*; do
  [ -d "${path}" ] || continue # if not a directory, skip
  dirname="$(basename "${path}")"

  echo "Entering folder ${dirname}"
  echo "Do you want to continue with this technology? [y/N]"
  read -r input

  if [[ ${input} =~ ^[Yy]$ ]]; then

    cd "${path}" || exit
    sanity_check=true

    ## Download data ##

    if [ ! -f "alignment.bam" ]; then
      filename='alignment.txt'
      while read -r url; do
        wget "${url}" -O "alignment.bam"
      done < ${filename}
    fi

    if [ ! -f "alignment.bam.bai" ] && [ -f "alignment_index.txt" ] ; then
      {
        filename='alignment_index.txt'
        while read -r url; do
          wget "${url}" -O "alignment.bam.bai"
        done < ${filename}
      } || echo "No alignment index file found"

    fi

    if [ ! -f "variants.vcf" ]; then
      filename='variants.txt'
      while read -r url; do
        wget "${url}" -O "variants.vcf"
      done < ${filename}
    fi

    ## Build variants ##

    if [ ! -f "${individual}"_chr1.var ] && [ ! -f "${individual}"_1.var ]; then
      {
        echo "Building variants ..."
        python2 ${script_folder}get.variants_mod.py "${individual}" ./variants.vcf
        echo "Variants builded"
        echo
      } || { sanity_check=false; }
    fi

    # Some variants are saved as 'individual'_chr_*, while others such as 'individual'_*

    individual_chr="${individual}"_chr1.var

    if [ -f "${individual}"_1.var ] ; then

      individual_chr="${individual}"_1.var

      elif [ ! -f "${individual}"_chr1.var  ]; then
      sanity_check=false;

    fi

    ## Build sfi file ##

    if [ "${sanity_check}" = true ] && [ ! -f "${individual}.sfi" ]; then
      #
      # get.sfi_mod uses pysam to fetch content of bam alignment file. If something goes wrong, fallback
      # on "old" script get.sfi.py in pipeline with samtools view

      {
        {
          echo "Building sfi file using pysam ..."
          python ${script_folder}get.sfi_mod.py "${individual_chr}" alignment.bam > "${individual}".sfi
          echo "${individual}.sfi created"
          echo
        } ||  # catch
        {
          echo "Fallback! Building sfi file using samtools ..."
          samtools view alignment.bam | python2 ${script_folder}get.sfi.py "${individual_chr}" > "${individual}".sfi
          echo "${individual}.sfi created"
          echo
        } || { sanity_check=false; }
      }
    fi

    ## Build matrix and its transpose ##

    if [ "${sanity_check}" = true ] && [ ! -f "${individual}.matrix" ] ; then
      {
        echo "Building sfi matrix ..."
        python2 ${script_folder}get.matrix.py "${individual_chr}" "${individual}".sfi > "${individual}".matrix
        echo "${individual}.matrix created"
        echo
      } ||  { sanity_check=false; }
    fi

    if [ "${sanity_check}" = true ] && [ ! -f "${individual}.transpose" ] ; then
      {
        echo "Building sfi matrix transpose ..."
        python2 ${script_folder}get.transpos.py "${individual_chr}" "${individual}".matrix > "${individual}".transpos
        echo "${individual}.transpos created"
      } || { sanity_check=false; }
    fi

    ## Build variants stats ##

    if [ "${sanity_check}" = true ] && [ ! -f "${individual}.stats" ] ; then
      {
        echo "Building individual stats ..."
        python2 ${script_folder}get.var_stats_mod.py "${individual}".transpos > "${individual}".stats
        echo "${individual}.stats created"
        echo
      } || { sanity_check=false; }
    fi

    ## Delete some useless file ##

    if [ "${sanity_check}" = true ] ; then
      echo "All fine. Can remove useless data."
      echo "Do you want to clear bam, sfi and other chr data? [y/N]"
      read -r input

      if [[ ${input} =~ ^[Yy]$ ]]; then {
          rm -rf alignment.bam alignment.bam.bai variants.vcf "${individual}_*" "${individual}.sfi" "${individual}.matrix" "${individual}.transpos"
        }
      fi
    else
      echo "Something went wrong."
      echo
    fi

    echo "Computation finished for ${dirname}"
    echo
  fi
done
#!/usr/bin/env python3
# - * -coding: utf - 8 - * -

import logging
import os
import sys
import argparse
import progressbar
import re
import pprint

import pysam


def check_observed(observed_couple):
    # observed are in the form A/T, T/-, -/C
    regexpr_snp = re.compile(r'^(\-|\w+)\/(\-|\w+)')
    couple = regexpr_snp.search(observed_couple)
    return couple.group(1), couple.group(2)


def file_len(file):
    # return the length (number of rows) of a file
    with open(file) as f:
        for i, l in enumerate(f):
            pass
    return i


def snpsInRead(alignment_file, annotation_file, file_delimiter, skip_first_line):

    in_sam = pysam.AlignmentFile(alignment_file, 'rb')

    # return the number of SNPs in the file
    # snp_total_number = file_len(annotation_file)

    total_fetch_alignment = in_sam.fetch()
    total_read_number = in_sam.count()

    read_count = 0
    read_bar = progressbar.ProgressBar(maxval=total_read_number).start()
    snp_read_list = {}
    for read in total_fetch_alignment:
        read_count += 1
        read_bar.update(read_count)

        read_qname = read.qname
        read_start = read.reference_start
        read_end = read.reference_end
        logging.debug('read {0}, start: {1}, end: {2}'.format(
            read_qname, read_start, read_end))

        ann_file = open(annotation_file)

        # skip first line (comment if the file starts with the header)
        if skip_first_line:
            ann_file.readline()

        # Header info - first line
        header = ann_file.readline().split(file_delimiter)

        # the fields to read
        snp_header_index_chromosome = header.index('chrom')
        snp_header_index_start = header.index('chromStart')
        snp_header_index_stop = header.index('chromEnd')
        snp_header_index_observed = header.index('observed')
        snp_read_list[read_qname] = []
        # snp_bar = progressbar.ProgressBar(maxval=snp_total_number).start()
        # snp_count = 0
        for line in ann_file:
            # snp_count += 1
            # snp_bar.update(snp_count)

            splitted_line = line.split(file_delimiter)

            snp_chrom = splitted_line[snp_header_index_chromosome]
            snp_start = int(splitted_line[snp_header_index_start])
            snp_stop = int(splitted_line[snp_header_index_stop])
            first_observed, second_observed = check_observed(
                splitted_line[snp_header_index_observed])

            # Check every read aligned with the SNPs
            for pileupcolumn in in_sam.pileup(snp_chrom, snp_start, snp_stop):
                for pileupread in pileupcolumn.pileups:
                    read_snp_name = pileupread.alignment.query_name
                    logging.debug('Found read: {0}'.format(read_snp_name))
                    if read_snp_name == read_qname:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            read_snp_position = pileupread.query_position
                            base_value = pileupread.alignment.query_sequence[
                                pileupread.query_position]
                            if base_value == first_observed:
                                allele = '0'
                            elif base_value == second_observed:
                                allele = '1'
                            else:
                                # skip the SNP if not in first or second allele
                                continue
                            base_quality = pileupread.alignment.query_alignment_qualities[
                                read_snp_position]
                            # if read_snp_name == read_qname and
                            # read_snp_position >= read_start and
                            # read_snp_position <= read_end:

                            snp_read_list[read_qname].append(
                                [read_snp_position, base_value, allele, base_quality])

        ann_file.close()
        snp_read_list[read_qname].append(read.mapping_quality)
        # snp_bar.finish()
    read_bar.finish()
    return snp_read_list


def convertBamToWif(snp_read_list, out_dir):

    output_file = open(out_dir + 'alignment.wif', 'w')

    for read, snp_list in snp_read_list.items():
        if len(snp_list) <= 1:
            continue
        else:
            snp_mapping_quality = snp_list[-1]
            for snp in snp_list[:-1]:
                output_file.write('{0} {1} {2} {3} : '.format(
                    snp[0], snp[1], snp[2], snp[3]))
            output_file.write('# {0} : NA\n'.format(snp_mapping_quality))


def main():

    parser = argparse.ArgumentParser(
        description='Convert an alignment file .bam in a .wif Hapcol compatible file')
    parser.add_argument('-b',
                        action='store',
                        dest='fileBAM',
                        help='Alignment file in BAM format.',
                        required=True)
    parser.add_argument('-a',
                        action='store',
                        dest='annotationFile',
                        help='SNPs annotation file.',
                        required=True)
    parser.add_argument('-o',
                        action='store',
                        dest='output_dir',
                        help='Output (root) directory.',
                        required=True)
    parser.add_argument('-skip-first-line',
                        action='store',
                        dest='skip_first_line',
                        help='set this flag if the first line IS NOT the header')
    parser.add_argument('-file-delimiter',
                        action='store',
                        dest='file_delim',
                        help='set the file delimiter for the annotation file. Default: \\t')
    parser.add_argument('-v',
                        help='increase output verbosity',
                        action='count')

    args = parser.parse_args()

    logging.info('#### STEP 1 - PROCESSING INPUT ARGUMENTS ####')

    in_bam_file = args.fileBAM
    out_root_dir = args.output_dir

    if args.v == 0:
        log_level = logging.INFO
    else:
        log_level = logging.DEBUG

    logging.basicConfig(level=log_level,
                        format='%(levelname)-8s [%(asctime)s]  %(message)s',
                        datefmt='%y%m%d %H%M%S')

    if not os.path.exists(out_root_dir):
        logging.error('Output dir not found.')
        sys.exit(1)

    logging.info('BamToWif: Program Started')

    out_dir = out_root_dir + '/'

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    if args.file_delim is None:
        ann_file_delimiter = '\t'
    else:
        ann_file_delimiter = args.file_delimiter

    logging.info('#### STEP 2 - RETREIVING SNPs FOR READS ####')

    snps_list = snpsInRead(in_bam_file, args.annotationFile,
                           ann_file_delimiter, args.skip_first_line)

    logging.info('#### STEP 3 - CONVERTING BAM TO WIF ####')

    convertBamToWif(snps_list, out_dir)

    logging.info('BamToWif: Program Completed')

if __name__ == '__main__':
    main()

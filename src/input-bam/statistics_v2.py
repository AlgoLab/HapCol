#!/usr/bin/env python3
# - * -coding: utf - 8 - * -

import argparse
import logging
import os
import re
import sys

import progressbar
import pysam
import time


def generate_statistics(snp_map, alignment_file):
    logging.debug('Started method generate_statistics')
    start_time = time.time()

    in_sam = pysam.AlignmentFile(alignment_file, 'rb')

    snp_total_number = len(snp_map)

    snp_count = 0
    prgbar = progressbar.ProgressBar(maxval=snp_total_number).start()
    count_first_allele = 0
    count_second_allele = 0
    count_other_allele = 0

    read_count = 0
    for snp_name, snp_values in snp_map.items():
        snp_count += 1
        prgbar.update(snp_count)
        logging.debug('Loaded snp: {0}'.format(snp_name))
        snp_chrom = snp_values[0]
        snp_start = snp_values[1]
        snp_end = snp_values[2]
        first_observed = snp_values[3]
        second_observed = snp_values[4]
        logging.debug(
            'Reading chromosome: {0}, Starting: {1} - Ending{2}'.format(snp_chrom, snp_start, snp_end))
        fetch_aln = in_sam.fetch(snp_chrom, snp_start, snp_end)
        tot_fetch_aln = in_sam.count(snp_chrom, snp_start, snp_end)
        logging.info(
            'Number of reads aligned with snp: {0}'.format(tot_fetch_aln))

        if tot_fetch_aln == 0 or tot_fetch_aln < 10:
            logging.debug('SNP without at least ten items will be discarded.')
            continue

        for read in fetch_aln:
            read_count += 1
            read_name = read.query_name
            read_start = read.reference_start
            read_end = read.reference_end
            read_sequence = read.query_sequence
            logging.debug('Loaded read: {0}'.format(read_name))

            if snp_start >= read_start and snp_end <= read_end:
                logging.debug('In the read: {0}, in the SNP (first): {1}, in the SNP (second): {2}'.format(
                    read_sequence[snp_start:snp_end], first_observed, second_observed))
                if read_sequence[snp_start:snp_end] == first_observed:
                    count_first_allele += 1
                elif read_sequence[snp_start:snp_end] == second_observed:
                    count_second_allele += 1
                elif not read_sequence[snp_start:snp_end].strip():
                    logging.debug('No valid read found.')
                else:
                    count_other_allele += 1
    prgbar.finish()

    logging.info('Total number of SNPs: {0}'.format(snp_total_number))

    logging.info('Total number of reads: {0}'.format(read_count))

    logging.info('Total number of read aligned with the first allele: {0} ({1}%)'.format(
        count_first_allele, (float(count_first_allele) / read_count) * 100))

    logging.info('Total number of read aligned with the second allele: {0} ({1}%)'.format(
        count_second_allele, (float(count_second_allele) / read_count) * 100))

    logging.info('Total number of read aligned with other allele: {0} ({1}%)'.format(
        count_other_allele, (float(count_other_allele) / read_count) * 100))

    logging.debug('Finished method generate_statistics in {0} seconds'.format(
        time.time() - start_time))
    logging.info('Generation of statistics: Completed')
    in_sam.close()


def main():

    parser = argparse.ArgumentParser(
        description='Generate statistics using alignment BAM file and SNPs annotation file')
    parser.add_argument('-b',
                        action='store',
                        dest='fileBAM',
                        help='Alignment file in BAM format.',
                        required=True)
    parser.add_argument('-a',
                        action='store',
                        dest='annotationFile',
                        help='SNPs location annotation file.',
                        required=True)
    parser.add_argument('-o',
                        action='store',
                        dest='output_dir',
                        help='Output (root) directory.',
                        required=True)
    parser.add_argument('-skip-first-line',
                        action='store_true',
                        dest='skip_first_line',
                        help='set this flag if the first line IS NOT the header')
    parser.add_argument('-file-delimiter',
                        action='store',
                        dest='file_delim',
                        help='set the file delimiter for the annotation file. Default: \\t')
    parser.add_argument('-v',
                        help='increase output verbosity',
                        dest='verbosity',
                        action='store_true')

    args = parser.parse_args()

    in_bam_file = args.fileBAM
    out_root_dir = args.output_dir

    if args.verbosity:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    if not os.path.exists(out_root_dir):
        logging.error('Output dir not found.')
        sys.exit(1)

    out_dir = out_root_dir + '/'

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    prepare_loggers(log_level)

    logging.info('Program Started')

    if args.file_delim is None:
        ann_file_delimiter = '\t'
    else:
        ann_file_delimiter = args.file_delimiter

    logging.info('#### STEP 1 - Loading SNPs ####')

    snp_map = reading_snp(args.annotationFile,
                          ann_file_delimiter, args.skip_first_line)

    logging.info('#### STEP 2 - Generating Statistics ####')

    generate_statistics(snp_map, in_bam_file)

    logging.info('Program Completed')


# utility

def check_observed(observed_couple):
    # observed are in the form A/T, T/-, -/C
    regexpr_snp = re.compile(r'^(\-|\w+)\/(\-|\w+)')
    couple = regexpr_snp.search(observed_couple)
    return couple.group(1), couple.group(2)


def prepare_loggers(log_level):
    logging.basicConfig(level=log_level,
                        format='%(levelname)-8s [%(asctime)s]  %(message)s',
                        datefmt='%y%m%d %H:%M:%S')


def reading_snp(annotation_file, file_delimiter, flag_first_line):

    logging.debug('Started method reading_snp')
    start_time = time.time()
    ann_file = open(annotation_file)

    # skip first line (comment if the file starts with the header)
    if flag_first_line:
        ann_file.readline()

    # read header info
    header = ann_file.readline().split(file_delimiter)

    # the fields to read
    snp_header_index_chrom = header.index('chrom')
    snp_header_index_start = header.index('chromStart')
    snp_header_index_stop = header.index('chromEnd')
    snp_header_index_name = header.index('name')
    snp_header_observed = header.index('observed')
    snp_map = {}
    for line in ann_file:
        splitted_line = line.split(file_delimiter)

        snp_chrom = splitted_line[snp_header_index_chrom].split('chr', 1)[1]
        snp_start = int(splitted_line[snp_header_index_start])
        snp_stop = int(splitted_line[snp_header_index_stop])
        snp_name = splitted_line[snp_header_index_name]
        first_observed, second_observed = check_observed(
            splitted_line[snp_header_observed])
        logging.debug('{0},{1},{2},{3}'.format(
            snp_chrom, snp_start, snp_stop, snp_name))
        snp_map[snp_name] = [snp_chrom, snp_start, snp_stop,
                             first_observed, second_observed]

    ann_file.close()

    logging.debug('Finished method reading_snp in: {0} seconds'.format(
        time.time() - start_time))

    return snp_map


if __name__ == '__main__':
    main()

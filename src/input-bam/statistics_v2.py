#!/usr/bin/env python3
# - * -coding: utf - 8 - * -

import argparse
import logging
import os
import re
import sys
from collections import Counter
import time

import progressbar
import pysam


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
    parser.add_argument('-epsilon',
                        action='store',
                        dest='threshold',
                        help='set the threshold for Homozygous or Heterozygous. Default: 80')
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

    threshold = 80
    if args.threshold is not None:
        threshold = int(args.threshold)

    logging.info('#### STEP 1 - Loading SNPs ####')

    snp_map = reading_snp(args.annotationFile,
                          ann_file_delimiter, args.skip_first_line)

    logging.info('#### STEP 2 - Generating Statistics ####')

    classify_snp(in_bam_file, snp_map, threshold, out_dir)

    logging.info('Program Completed')


# utility
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
    snp_map = {}
    for line in ann_file:

        splitted_line = line.split(file_delimiter)

        snp_chrom = splitted_line[snp_header_index_chrom].split('chr', 1)[1]
        snp_start = int(splitted_line[snp_header_index_start])
        snp_stop = int(splitted_line[snp_header_index_stop])
        snp_name = splitted_line[snp_header_index_name]

        snp_map[snp_name] = [snp_chrom, snp_start, snp_stop]

    ann_file.close()

    logging.debug('Finished method reading_snp in: {0} seconds'.format(
        time.time() - start_time))
    return snp_map


def classify_snp(alignment_file, snp_map, threshold, out_dir):
    logging.debug('Started method classifySNP')
    start_time = time.time()
    count_homo = 0
    count_hete = 0
    snp_total_number = len(snp_map)
    statistics_file = open(out_dir + 'statistics2.txt', 'w')

    in_sam = pysam.AlignmentFile(alignment_file, 'rb')
    prgbar = progressbar.ProgressBar(maxval=snp_total_number).start()
    cnt = 0

    for snp_name, snp_values in snp_map.items():
        statistics_file.write('Reading SNP: {0}\n'.format(snp_name))
        cnt += 1
        prgbar.update(cnt)
        snp_chrom = snp_values[0]
        snp_start = snp_values[1]
        snp_end = snp_values[2]

        logging.debug(
            'SNP start: {0} - SNP end:{1}'.format(snp_start, snp_end))

        fetch_aln = in_sam.fetch(snp_chrom, snp_start, snp_end)
        count_aln = in_sam.count(snp_chrom, snp_start, snp_end)

        alleles = Counter()

        for read in fetch_aln:
            read_name = read.query_name
            logging.debug('Reading read: {0}'.format(read_name))
            read_sequence = read.query_sequence
            read_base = read_sequence[snp_start:snp_end]
            if not read_base:
                logging.debug('Empty read - indel')
                continue
            logging.debug('Read base: {0}'.format(read_base))
            alleles[read_base] += 1

        logging.debug(alleles)
        statistics_file.write('Allele counter for the SNP:\n')
        for key, value in alleles.items():
            statistics_file.write('{0} : {1}\n'.format(key, value))

        if alleles:
            first_allele = alleles.most_common(2)[0]
            logging.debug(first_allele)
            second_allele = alleles.most_common(2)[1]
            logging.debug(second_allele)

            if (first_allele[1] / count_aln) * 100 >= threshold:
                statistics_file.write('Snp is Homozygous\n')
                logging.debug('The snp {0} is: Homozygous'.format(snp_name))
                count_homo += 1
            else:
                statistics_file.write('Snp is Heterozygous\n')
                logging.debug('The snp {0} is: Heterozygous'.format(snp_name))
                count_hete += 1

        logging.debug(
            'Number of read aligned on SNP: {0}'.format(count_aln))
        statistics_file.write(
            'Number of read aligned on SNP: {0}\n'.format(count_aln))

        out = {}
        if count_aln > 0:
            for key in alleles:
                out[key] = (float(alleles[key]) / count_aln) * 100

        logging.debug(out)
        statistics_file.write('Base Frequences for the SNP:\n')
        for key, value in out.items():
            statistics_file.write('{0} : {1}\n'.format(key, value))

    prgbar.finish()
    in_sam.close()

    logging.debug(
        'Total number of SNPs: {0}'.format(snp_total_number))

    statistics_file.write(
        'Total number of SNPs: {0}\n'.format(snp_total_number))

    logging.debug(
        'Number of SNPs: Homozygous {0} ({1}%)'.format(
            count_homo, ((float(count_homo) / snp_total_number) * 100)))

    statistics_file.write('Number of SNPs: Homozygous {0} ({1}%)\n'.format(
        count_homo, ((float(count_homo) / snp_total_number) * 100)))

    logging.debug(
        'Number of SNPs: Heterozygous {0} ({1}%)'.format(
            count_hete, ((float(count_hete) / snp_total_number) * 100)))

    statistics_file.write('Number of SNPs: Heterozygous {0} ({1}%)\n'.format(
        count_hete, ((float(count_hete) / snp_total_number) * 100)))

    logging.debug('Finished method classifySNP in: {0} seconds'.format(
        time.time() - start_time))

    statistics_file.close()


if __name__ == '__main__':
    main()

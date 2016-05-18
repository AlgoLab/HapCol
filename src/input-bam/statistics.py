#!/usr/bin/env python3
# - * -coding: utf - 8 - * -

import argparse
import logging
import os
import re
import sys

import progressbar
import pysam


def generate_statistics(annotation_file, file_delimiter, region, alignment_file, out_dir, skip_first_line):
    in_sam = pysam.AlignmentFile(alignment_file, 'rb')

    # return the number of SNPs in the file
    snp_total_number = file_len(annotation_file)
    logging.info('Number of Total SNPs: {0}'.format(snp_total_number))
    ann_file = open(annotation_file)

    # skip first line (comment if the file starts with the header)
    if skip_first_line:
        ann_file.readline()

    # read header info
    header = ann_file.readline().split(file_delimiter)

    # the fields to read
    snp_header_index_start = header.index('chromStart')
    snp_header_index_stop = header.index('chromEnd')
    snp_header_index_name = header.index('name')
    snp_header_observed = header.index('observed')

    statistics_file = open(out_dir + 'statistics.txt', 'w')

    # read annotation file
    snp_valid = 0  # SNPs with almost 1 read aligned
    snp_not_valid = 0  # SNPs with no read aligned
    snp_count = 0
    read_with_first_allele = 0
    read_with_second_allele = 0
    read_with_different_allele = 0
    prgbar = progressbar.ProgressBar(maxval=snp_total_number).start()
    for line in ann_file:
        snp_count += 1
        prgbar.update(snp_count)

        splitted_line = line.split(file_delimiter)

        snp_start = int(splitted_line[snp_header_index_start])
        snp_stop = int(splitted_line[snp_header_index_stop])
        snp_name = splitted_line[snp_header_index_name]
        first_observed, second_observed = check_observed(
            splitted_line[snp_header_observed])

        frequences_observed = {'first': 0, 'second': 0, 'others': 0}

        logging.info('Loaded snp: {0}'.format(snp_name))
        statistics_file.write('\nLoaded snp: {0}'.format(snp_name))
        frequences = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
        count_read = 0
        for pileupcolumn in in_sam.pileup(region[0], snp_start, snp_stop):

            logging.debug('Coverage at base {0} = {1}'.format(
                pileupcolumn.pos, pileupcolumn.n))
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    count_read += 1
                    base_value = pileupread.alignment.query_sequence[
                        pileupread.query_position]
                    logging.debug('Base in read = {0}, First observed = {1}, Second observed = {2}'.format(
                        base_value, first_observed, second_observed))
                    frequences[base_value] += 1
                    if base_value == first_observed:
                        frequences_observed['first'] += 1
                        read_with_first_allele += 1
                    elif base_value == second_observed:
                        frequences_observed['second'] += 1
                        read_with_second_allele += 1
                    else:
                        frequences_observed['others'] += 1
                        read_with_different_allele += 1
        logging.debug(frequences)

        if count_read != 0:
            snp_valid += 1
            # number of reads aligned with the SNP
            logging.info(
                'Total number of read aligned with SNP: {0}'.format(count_read))
            statistics_file.write(
                '\nTotal number of read aligned with SNP: {0}'.format(count_read))

            # number of reads aligned with the SNP in the first allele
            logging.info('Number of read aligned with the first allele: {0}'.format(
                frequences_observed['first']))
            statistics_file.write('\nNumber of read aligned with the first allele: {0}'.format(
                frequences_observed['first']))

            # reads aligned with the SNP in the first allele
            logging.info('Percentual of read aligned with the first allele: {0}'.format(float(
                frequences_observed['first']) / count_read * 100))
            statistics_file.write('\nPercentual of read aligned with the first allele: {0}'.format(float(
                frequences_observed['first']) / count_read * 100))

            # number of reads aligned with the SNP in the second allele
            logging.info('Number of read aligned with the second allele: {0}'.format(
                frequences_observed['second']))
            statistics_file.write('\nNumber of read aligned with the second allele: {0}'.format(
                frequences_observed['second']))

            # reads aligned with the SNP in the second allele
            logging.info('Percentual of read aligned with the second allele: {0}'.format(float(
                frequences_observed['second']) / count_read * 100))
            statistics_file.write('\nPercentual of read aligned with the second allele: {0}'.format(float(
                frequences_observed['second']) / count_read * 100))

            # reads aligned with the SNP - different allele
            logging.info('Number of read aligned with different allele: {0}'.format(
                frequences_observed['others']))
            statistics_file.info('\nNumber of read aligned with different allele: {0}'.format(
                frequences_observed['others']))

            logging.info('Percentual of read aligned with different allele: {0}'.format(
                float(frequences_observed['others']) / count_read) * 100)
            statistics_file.write('\nPercentual of read aligned with different allele: {0}'.format(
                float(frequences_observed['others']) / count_read) * 100)

        else:
            snp_not_valid += 1
            logging.info(
                'No reads aligned with this snp: {0}'.format(snp_name))
            statistics_file.write(
                '\nNo reads aligned with this snp: {0}'.format(snp_name))

    prgbar.finish()
    logging.info(
        '#### TOTAL RESULTS ####')
    statistics_file.write(
        '\n#### TOTAL RESULTS ####')

    logging.info('Processed {0} SNPs'.format(snp_count))
    statistics_file.write('\nProcessed {0} SNPs'.format(snp_count))

    logging.info('SNPs with almost one read aligned (%): {0}'.format(
        float(snp_valid) / snp_count) * 100)
    statistics_file.write('\nSNPs with almost one read aligned (%): {0}'.format(
        float(snp_valid) / snp_count) * 100)
    logging.info('SNPs with no read aligned (%): {0}'.format(
        float(snp_not_valid) / snp_count) * 100)
    statistics_file.write('\nSNPs with no read aligned (%): {0}'.format(
        float(snp_not_valid) / snp_count) * 100)

    logging.info('Total number of read aligned with the first allele: {0}'.format(
        read_with_first_allele))
    statistics_file.write(
        '\nTotal number of read aligned with the first allele: {0}'.format(read_with_first_allele))

    logging.info('Total number of read aligned with the second allele: {0}'.format(
        read_with_second_allele))
    statistics_file.write(
        '\nTotal number of read aligned with the second allele: {0}'.format(read_with_second_allele))
    logging.info('Total number of read aligned with different allele: {0}'.format(
        read_with_different_allele))
    statistics_file.write(
        '\nTotal number of read aligned with different allele: {0}'.format(read_with_different_allele))

    logging.info('Generation of statistics: Completed')
    statistics_file.close()
    ann_file.close()
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
    parser.add_argument('-r',
                        action='store',
                        dest='region',
                        help='Region to be examined in the form: chr:start-stop. Ex. 1:100-200.')
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

    in_bam_file = args.fileBAM
    in_region = args.region
    out_root_dir = args.output_dir

    if args.v is None:
        log_level = logging.INFO
    else:
        log_level = logging.DEBUG

    logging.basicConfig(level=log_level,
                        format='%(levelname)-8s [%(asctime)s]  %(message)s',
                        datefmt='%y%m%d %H%M%S')

    if not os.path.exists(out_root_dir):
        logging.error('Output dir not found.')
        sys.exit(1)

    logging.info('StatisticsGenerator: Program Started')

    region = check_region(in_region)

    out_dir = out_root_dir + '/'

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    if args.file_delim is None:
        ann_file_delimiter = '\t'
    else:
        ann_file_delimiter = args.file_delimiter

    generate_statistics(args.annotationFile,
                        ann_file_delimiter, region, in_bam_file, out_dir, args.skip_first_line)

    logging.info('StatisticsGenerator: Program Completed')


# utility

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


def check_region(input_region):
    # if no region is set, is used all the 1st chromosome (for test purpose)
    if input_region is not None:
        regexp_reg = re.compile(r'^(\w+):(\d+)-(\d+)')
        region = re.search(regexp_reg, input_region)
        if region is not None:
            in_chr = region.group(1)
            in_start = int(region.group(2))
            in_stop = int(region.group(3))
            if in_start > in_stop:
                logging.error('Wrong input region.')
                sys.exit(1)
            else:
                return in_chr, in_start, in_stop
    else:
        return '1'


if __name__ == '__main__':
    main()

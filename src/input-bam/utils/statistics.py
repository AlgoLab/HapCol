#!/usr/bin/env python3
# - * -coding: utf - 8 - * -

import argparse
import logging
import os
import time

from tqdm import tqdm
import sys


def main():
    """
    It takes in input a variant statistics file ( in output from get.var_stats.py ) and returns how many SNPs
    heterozygous and homozygous are present using a threshold value
    :return: Info in stdout
    """
    parser = argparse.ArgumentParser(description='SNPs classification using a threshold value. Print stats to stdout')
    parser.add_argument('-i', action='store', dest='input_file', help='Variant statistics file', required=True)
    parser.add_argument('-e', action='store', dest='threshold', help='Threshold for classification. Default: 80',
                        type=int)
    parser.add_argument('-d', action='store', dest='file_del', help='Input file delimiter. Default: \\t')
    parser.add_argument('-o', action='store', dest='outputDir',
                        help='Output (root) directory. Default: current directory')
    parser.add_argument('-v', help='increase output verbosity', dest='verbosity', action='store_true')

    args = parser.parse_args()

    input_file = args.input_file

    if args.verbosity:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    prepare_loggers(log_level)

    logging.info('Program Started')

    if args.outputDir and not os.path.exists(args.outputDir):
        logging.error('Output dir not found.')
        sys.exit(1)

    if args.outputDir:
        out_dir = args.outputDir + '/'
    else:
        out_dir = os.getcwd() + '/'

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    if args.file_del is None:
        file_delimiter = '\t'
    else:
        file_delimiter = args.file_delim

    threshold = 80
    if args.threshold is not None:
        threshold = int(args.threshold)

    logging.info('#### Generating Statistics ####')

    classify_snp(input_file, file_delimiter, threshold, out_dir)

    logging.info('Program Completed')


# utility
def prepare_loggers(log_level):
    logging.basicConfig(level=log_level,
                        format='%(levelname)-8s [%(asctime)s]  %(message)s',
                        datefmt='%y%m%d %H:%M:%S')


def get_file_length(file):
    """
    Return the number of row ( header excluded )
    :param file:
    :return:
    """
    length = 0
    if file is not None:
        try:
            f = open(file, 'r')
            for line in f:
                if line.startswith('#'):
                    continue
                else:
                    length += 1
            f.close()
        except Exception as e:
            logging.error("Error during opening file. {0}".format(e))
    return length


def check_percentage(base_count, snp_number, threshold):
    return (float(base_count) / snp_number) * 100 >= threshold


def classify_snp(input_file, file_delimiter, threshold, out_dir):
    """
    :param input_file: individual.stats - from get.var_stats.py script
    :param file_delimiter: default \t
    :param threshold: default 80% - depends on error rate for technology
    :param out_dir:
    :return:
    """
    logging.debug('Started method classifySNP')
    start_time = time.time()
    snp_total_number = get_file_length(input_file)

    stats_file = open(input_file, 'r')

    count_homo = 0
    count_het = 0

    out_file = open(out_dir + 'homozygous.pos', 'w')
    for line in tqdm(stats_file, total=snp_total_number):
        if line.strip().startswith("#") or not line.strip():
            continue

        row = line.split(file_delimiter)

        homozygous_flag = False
        snp_position = row[0]
        snp_ref = row[1]
        snp_alt = row[2]
        snp__a = float(row[3])
        snp__c = float(row[5])
        snp__g = float(row[7])
        snp__t = float(row[9])
        snp_total_reads = float(row[12])

        if snp_total_reads > 0:
            logging.debug('Reading SNP: {} \n'.format(snp_position))
            if ((snp_ref == 'A' or snp_alt == 'A') and check_percentage(snp__a, snp_total_reads, threshold)) \
                    or check_percentage(snp__a, snp_total_reads, threshold):
                homozygous_flag = True
            elif ((snp_ref == 'C' or snp_alt == 'C') and check_percentage(snp__c, snp_total_reads, threshold)) \
                    or check_percentage(snp__c, snp_total_reads, threshold):
                homozygous_flag = True
            elif ((snp_ref == 'G' or snp_alt == 'G') and check_percentage(snp__g, snp_total_reads, threshold)) \
                    or check_percentage(snp__g, snp_total_reads, threshold):
                homozygous_flag = True
            elif ((snp_ref == 'T' or snp_alt == 'T') and check_percentage(snp__t, snp_total_reads, threshold)) \
                    or check_percentage(snp__t, snp_total_reads, threshold):
                homozygous_flag = True

            if homozygous_flag:
                count_homo += 1
                out_file.write("{}\n".format(snp_position))
                logging.debug("SNP {} is homozygous".format(snp_position))
            else:
                count_het += 1
                logging.debug("SNP {} is heterozygous".format(snp_position))

    stats_file.close()

    logging.info("Number or total SNPs: {} \n".format(snp_total_number))
    logging.info(
        "Number of Homozygous SNPs: {}, {}% \n".format(count_homo, (float(count_homo) / snp_total_number) * 100))
    logging.info(
        "Number of Heterozygous SNPs: {0}, {1}% \n".format(count_het, (float(count_het) / snp_total_number) * 100))

    logging.debug('Finished method classifySNP in: {} seconds'.format(
        time.time() - start_time))


if __name__ == '__main__':
    main()

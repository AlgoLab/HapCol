#!/usr/bin/env python3
# - * -coding: utf - 8 - * -

import argparse
import csv
import logging
import os
import sys

from pysam import VariantFile


# utility
def prepare_loggers(log_level):
    logging.basicConfig(level=log_level,
                        format='%(levelname)-8s [%(asctime)s]  %(message)s',
                        datefmt='%y%m%d %H:%M:%S')


def extract_accuracy_range(bed_file, chromosome):
    """

    :param bed_file:
    :param chromosome: if None, all chromosomes are read
    :return: an array containing tuples [(a,b),(b+1,c), (d,e) ...]
    """
    out = []
    with open(bed_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if chromosome and row[0] == chromosome:
                out.append((int(row[1]), int(row[2])))
            elif not chromosome:
                out.append((int(row[1]), int(row[2])))
    return out


def check_accuracy(position, ranges):
    """

    :param position:
    :param ranges: an array of tuples containing ranges --> [ (a,b), (c,d), ... ]
    :return: true if position is within one of the ranges
    """
    for r in ranges:
        r0, r1 = r
        if int(position) < r0:
            return False
        if r0 <= int(position) <= r1:
            return True
    return False


def filter_vcf(vcf_file, bed_file, out_dir, chromosome=None):
    ranges = extract_accuracy_range(bed_file, chromosome)

    vcf_in = VariantFile(vcf_file)  # auto-detect input format
    if chromosome:
        filename = out_dir + "chr{}_high_accuracy_positions.vcf".format(chromosome)
    else:
        filename = out_dir + "high_accuracy_positions.vcf"

    vcf_out = VariantFile(filename, 'w', header=vcf_in.header)

    for row in vcf_in:
        if chromosome and chromosome == row.chrom:
            if check_accuracy(row.pos, ranges):
                vcf_out.write(row)
        else:
            if check_accuracy(row.pos, ranges):
                vcf_out.write(row)


def main():
    parser = argparse.ArgumentParser(description='Given a vcf file in input and a bed file containing high accuracy '
                                                 'ranges, return a vcf with only high accuracy positions')

    parser.add_argument('-i', action='store', dest='input_file', help='Vcf file', required=True)
    parser.add_argument('-b', action='store', dest='bed_file', help='bed file with high accuracy ranges', required=True)
    parser.add_argument('-c', action='store', dest='chromosome', help='The chromosome to analyze, otherwise all '
                                                                      'chromosomes')
    parser.add_argument('-o', action='store', dest='outputDir',
                        help='Output (root) directory. Default: current directory')
    parser.add_argument('-v', help='increase output verbosity', dest='verbosity', action='store_true')

    args = parser.parse_args()

    if args.verbosity:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    prepare_loggers(log_level)

    if not os.path.exists(args.input_file):
        logging.error('Input file not exists.')
        sys.exit(1)

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

    filter_vcf(args.input_file, args.bed_file, out_dir, args.chromosome)
    logging.info("Program Finished")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
# - * -coding: utf - 8 - * -

import argparse
import csv
import logging
import os
import sys

from pysam import VariantFile
from tqdm import tqdm


# utility

def prepare_loggers(log_level):
    logging.basicConfig(level=log_level,
                        format='%(levelname)-8s [%(asctime)s]  %(message)s',
                        datefmt='%y%m%d %H:%M:%S')


def count_positions(bed_file):
    out = 0
    with open(bed_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            splitted_line = line.split('\t')
            r0 = int(splitted_line[1])
            r1 = int(splitted_line[2])
            out += r1 - r0
    return out


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def extract_accuracy_range(bed_file, chromosome):
    """

    :param bed_file:
    :param chromosome: if None, all chromosomes are read
    :return: the set of all high-confidence positions included in the bed file
    """
    out = set()
    with open(bed_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if chromosome and row[0] == chromosome or not chromosome:
                [out.add(x) for x in range(int(row[1]), int(row[2]) + 1)]
    return out


def add_range(range_list, chromosome, idx):
    out = set()
    for i, el in enumerate(range_list):
        if idx and i == idx:
            if chromosome and el[0] == chromosome or not chromosome:
                [out.add(x) for x in range(int(el[1]), int(el[2]) + 1)]
    return out


def get_range(range_list, chromosome, idx):
    for i, el in enumerate(range_list):
        if idx and i == idx:
            if chromosome and el[0] == chromosome or not chromosome:
                return el[1], el[2]


def file_to_list(file):
    """

    :param file: a bed file containing Chromosome RangeStart RangeEnd
    :return: a List representing the file's data as tuples [ (chr, start, end), (chr, start, end), ...]
    """
    out = []
    with open(file, 'r') as f:
        for line in f.readlines():
            row = line.split('\t')
            t = (row[0], int(row[1]), int(row[2]))
            out.append(t)
    return out


def filter_vcf(vcf_file, bed_file, out_dir, chromosome=None, snps=None):
    vcf_in = VariantFile(vcf_file)  # auto-detect input format
    if chromosome and snps:
        filename = out_dir + "chr{}_high_accuracy_positions_{}.vcf".format(chromosome, snps)
    elif chromosome:
        filename = out_dir + "chr{}_high_accuracy_positions.vcf".format(chromosome)
    elif snps:
        filename = out_dir + "high_accuracy_positions_{}.vcf".format(snps)
    else:
        filename = out_dir + "high_accuracy_positions.vcf"

    vcf_out = VariantFile(filename, 'w', header=vcf_in.header)

    if not snps:
        # all high accuracy positions are inserted into a set
        ranges = extract_accuracy_range(bed_file, chromosome)
        # then, check for every SNP if it's included in one of the high accuracy range.
        for row in vcf_in:
            if chromosome and chromosome == row.chrom or not chromosome:
                if row.pos in ranges:
                    vcf_out.write(row)
    else:
        snps_found = 0

        # is easier to work on ranges list instead file lines
        range_list = file_to_list(bed_file)
        ranges_added = 0
        useless_range = 0

        if chromosome and snps:
            filename = out_dir + "chr{}_high_accuracy_filtered_{}.bed".format(chromosome, snps)
        else:
            filename = out_dir + "high_accuracy_filtered_{}.bed".format(chromosome)

        bed_out = open(filename, 'w')

        if chromosome and snps:
            filename = out_dir + "chr{}_discarded_ranges_{}.bed".format(chromosome, snps)
        else:
            filename = out_dir + "discarded_ranges_{}.bed".format(chromosome)

        bed_discarded = open(filename, 'w')

        for i in tqdm(range_list, total=len(range_list), desc='Number of ranges'):
            vcf_in = VariantFile(vcf_file)
            # the computation end if:
            # - number of snps found is greater than a threshold
            # - all ranges are read
            if snps_found >= snps:
                break

            # start from the middle of the list
            middle = len(range_list) // 2
            if middle >= 0 and range_list:
                # add range one by one and check how many snps from vcf are included
                ranges = add_range(range_list, chromosome, middle)
                ranges_added += 1

                for row in vcf_in:
                    if chromosome and chromosome == row.chrom or not chromosome:
                        if row.pos in ranges:
                            snps_found += 1
                            vcf_out.write(row)

                a, b = get_range(range_list, chromosome, middle)
                if snps_found == 0:
                    bed_discarded.write("{}\t{}\t{}\n".format(chromosome, a, b))
                    useless_range += 1
                else:
                    bed_out.write("{}\t{}\t{}\n".format(chromosome, a, b))
                # remove the middle element from the range list
                range_list.pop(middle)
            vcf_in.close()
        bed_out.close()

        logging.info("Found {} SNPs".format(snps_found))
        logging.info("Number of ranges discarded {}".format(useless_range))


def main():
    parser = argparse.ArgumentParser(description='Given a vcf file in input and a bed file containing high accuracy '
                                                 'ranges, return a vcf with only high accuracy positions')

    parser.add_argument('-i', action='store', dest='input_file', help='Vcf file', required=True)
    parser.add_argument('-b', action='store', dest='bed_file', help='bed file with high accuracy ranges', required=True)
    parser.add_argument('-c', action='store', dest='chromosome', help='The chromosome to analyze, otherwise all '
                                                                      'chromosomes')
    parser.add_argument('-n', action='store', dest='number_of_snps',
                        help='Limit the number of snps to insert in the vcf. Default= all snps are retrieved', type=int)
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

    l = file_len(args.bed_file)
    c = count_positions(args.bed_file)
    logging.info("Number of total ranges: {} - Number of total positions: {}".format(l, c))
    filter_vcf(args.input_file, args.bed_file, out_dir, args.chromosome, args.number_of_snps)
    logging.info("Program Finished")


if __name__ == '__main__':
    main()

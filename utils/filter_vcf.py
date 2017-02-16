#!/usr/bin/env python3
# - * -coding: utf - 8 - * -

import argparse
import csv
import logging
import os
import sys
from enum import Enum

from pysam import VariantFile
from tqdm import tqdm


# utility

def prepare_loggers(log_level):
    logging.basicConfig(level=log_level,
                        format='%(levelname)-8s [%(asctime)s]  %(message)s',
                        datefmt='%y%m%d %H:%M:%S')


class FileType(Enum):
    BED = "bed"
    VCF = "vcf"


def filename_builder(prefix, chromosome, snp, file_type, out_dir):
    out = out_dir + prefix
    if chromosome:
        out += "_chr{}".format(chromosome)
    if snp:
        out += "_snp{}".format(snp)
    if file_type == FileType.BED:
        out += ".{}".format(FileType.BED.value)
    elif file_type == FileType.VCF:
        out += ".{}".format(FileType.VCF.value)

    return out


def count_positions_in_ranges(bed_file):
    """
    :param bed_file: bed file containing high-confidence ranges
    :return: count of all positions inside ranges
    """
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
    i = -1
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


def add_range(range_list, chromosome, idx, ranges_set):
    for i, el in enumerate(range_list):
        if idx and i == idx:
            if chromosome and el[0] == chromosome or not chromosome:
                [ranges_set.add(x) for x in range(int(el[1]), int(el[2]) + 1)]
    return ranges_set


def get_range(range_list, chromosome, idx):
    for i, el in enumerate(range_list):
        if idx and i == idx:
            if chromosome and el[0] == chromosome or not chromosome:
                return el[1], el[2]


def file_to_list(fname):
    """
    :param fname: bed file containing Chromosome RangeStart RangeEnd
    :return: list containing each file's row as tuple [ (Chromosome, RangeStart, RangeEnd), ...]
    """
    out = []
    with open(fname, 'r') as f:
        for line in f.readlines():
            row = line.split('\t')
            t = (row[0], int(row[1]), int(row[2]))
            out.append(t)
    return out


def check_snv(ref, alts):
    if ref and alts:
        try:
            alt = alts[0]
            if len(ref) == len(alt) == 1:
                return True
        except IndexError:
            logging.exception("Alts tuple empty: {}".format(alts))
            return False
    return False


def filter_vcf(vcf_file, bed_file, out_dir, chromosome=None, snp_threshold=None, only_snps=True):
    vcf_out_filename = filename_builder("high_accuracy_snp", chromosome, snp_threshold, FileType.VCF, out_dir)

    vcf_header, vcf_length = extract_vcf_info(vcf_file)

    if not snp_threshold:
        vcf_in = VariantFile(vcf_file)  # auto-detect input format
        vcf_out = VariantFile(vcf_out_filename, 'w', header=vcf_header)
        # If no snp threshold is set, each position of bed file is read from the beginning.
        # All high confidence positions are inserted into a set
        ranges_set = extract_accuracy_range(bed_file, chromosome)
        # then, check for every vcf row check if row.pos ( position) is in high accuracy range.
        for row in tqdm(vcf_in, desc="Vcf rows", total=vcf_length):
            if chromosome and chromosome == row.chrom or not chromosome:
                if row.pos in ranges_set:
                    if (only_snps and check_snv(row.ref, row.alts)) or not only_snps:
                        vcf_out.write(row)
        vcf_out.close()
        vcf_in.close()

    else:
        snps_found = 0

        ranges_list = file_to_list(bed_file)
        ranges_set = set()

        vcf_out_filename = filename_builder("high_accuracy_snp", chromosome, snp_threshold, FileType.VCF, out_dir)
        bed_filename = filename_builder("high_accuracy_ranges", chromosome, snp_threshold, FileType.BED, out_dir)
        bed_discarded_filename = filename_builder("discarded_ranges", chromosome, snp_threshold, FileType.BED, out_dir)

        bed_out = open(bed_filename, 'w')
        bed_discarded = open(bed_discarded_filename, 'w')

        vcf_out = VariantFile(vcf_out_filename, 'w', header=vcf_header)

        for i in tqdm(ranges_list, total=len(ranges_list), desc='High accuracy ranges'):
            vcf_in = VariantFile(vcf_file)  # auto-detect input format
            # the computation end if:
            # - number of snps found is greater than a threshold
            # - all ranges are read
            if snps_found >= snp_threshold:
                break

            # start from the middle of the list
            middle_pos = len(ranges_list) // 2
            if middle_pos >= 0 and ranges_list:
                # add range one by one and check how many snps from vcf are included
                ranges_set = add_range(ranges_list, chromosome, middle_pos, ranges_set)

                for row in vcf_in:
                    if chromosome and chromosome == row.chrom or not chromosome:
                        if row.pos in ranges_set:
                            if (only_snps and check_snv(row.ref, row.alts)) or not only_snps:
                                snps_found += 1
                                vcf_out.write(row)

                a, b = get_range(ranges_list, chromosome, middle_pos)
                if snps_found == 0:
                    bed_discarded.write("{}\t{}\t{}\n".format(chromosome, a, b))
                else:
                    bed_out.write("{}\t{}\t{}\n".format(chromosome, a, b))
                # remove the middle_pos element from the range list
                ranges_list.pop(middle_pos)
            vcf_in.close()

        vcf_out.close()
        bed_discarded.close()
        bed_out.close()

        logging.info("Found {} SNPs".format(snps_found))
        d = file_len(bed_discarded_filename)
        logging.info("Number of ranges discarded {}".format(d))
        l = file_len(bed_filename)
        c = count_positions_in_ranges(bed_filename)
        logging.info("Number of founded ranges: {} - Number of positions: {}".format(l, c))


def extract_vcf_info(vcf_file):
    vcf_in = VariantFile(vcf_file)
    vcf_header = vcf_in.header
    vcf_rows = 0
    for row in vcf_in:
        vcf_rows += 1
    vcf_in.close()
    return vcf_header, vcf_rows


def main():
    parser = argparse.ArgumentParser(description='Given a vcf file in input and a bed file containing high accuracy '
                                                 'ranges, return a vcf with only high accuracy positions')

    parser.add_argument('-i', action='store', dest='input_file', help='Vcf file', required=True)
    parser.add_argument('-b', action='store', dest='bed_file', help='bed file with high accuracy ranges', required=True)
    parser.add_argument('-c', action='store', dest='chromosome', help='The chromosome to analyze, otherwise all '
                                                                      'chromosomes')
    parser.add_argument('-n', action='store', dest='number_of_snps',
                        help='Limit the number of snps to insert in the vcf. Default= all snps are retrieved', type=int)
    parser.add_argument('-s', action='store_true', dest='only_snp', help='Only SNPs. Default True', default=True)
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
    c = count_positions_in_ranges(args.bed_file)
    logging.info("Number of total ranges: {} - Number of total positions: {}".format(l, c))
    filter_vcf(args.input_file, args.bed_file, out_dir, args.chromosome, args.number_of_snps, args.only_snp)
    logging.info("Program Finished")


if __name__ == '__main__':
    main()

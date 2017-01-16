#!/usr/bin/env python3
# - * -coding: utf - 8 - * -

import argparse
import logging
import os
import sys

import pysam
from tqdm import tqdm


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def mean(list_of_quals):
    """
    Naive function for determining 'read quality'
    :param list_of_quals: list of numeric quality values
    :return: mean of the quality values
    """
    return float(sum(list_of_quals)) / len(list_of_quals)


def load_snps(alignment_file, variant_file, file_delimiter, chromosome, base_q=None, read_q=None):
    in_sam = pysam.AlignmentFile(alignment_file, 'rb')

    total_read_number = in_sam.count()
    logging.info('Number of fetched reads: {}'.format(total_read_number))

    out_map = {}
    total_variant_number = file_len(variant_file)
    var_file = open(variant_file)

    valid = ['A', 'T', 'C', 'G']

    x = 0
    for line in tqdm(var_file, total=total_variant_number, desc="Variant counter", miniters=1):

        splitted = line.split(file_delimiter)
        snp_position = int(splitted[0])
        snp_ref = splitted[1]
        snp_alt = splitted[2]

        for read in in_sam.fetch(chromosome, snp_position, snp_position + 1):

            if read_q and mean(read.query_qualities) < read_q:
                # read quality threshold
                continue

            if read.qname not in out_map:
                out_map[read.qname] = {'snp_list': [], 'mq': 0}

            try:
                base_value = read.seq[snp_position - read.reference_start]
            except IndexError:
                continue

            if base_value not in valid:
                continue

            if base_value == snp_ref:
                allele = 0
            elif base_value == snp_alt:
                allele = 1
            else:
                continue

            try:
                base_quality = read.query_qualities[snp_position - read.reference_start]
                if base_q and base_quality < base_q:
                    # base quality threshold
                    continue
            except IndexError:
                continue

            logging.debug("Found Snp: {} {} {} {}".format(snp_position, base_value, allele, base_quality))
            out_map[read.qname]['snp_list'].append((snp_position, base_value, allele, base_quality))
            out_map[read.qname]['mq'] = read.mapping_quality

    var_file.close()
    in_sam.close()
    return out_map


def write_wif(read_snp_map, out_dir):
    output_file = open(out_dir + 'alignment.wif', 'w')

    for read_name, snp_dict in read_snp_map.items():
        snp_list = snp_dict['snp_list']
        mq = snp_dict['mq']

        if len(snp_list) == 0:  # No SNP found for this read or read empty
            continue

        output_file.write("{}".format(read_name))
        for snp in snp_list:
            p = snp[0]  # position
            b = snp[1]  # read base
            a = snp[2]  # allele
            q = snp[3]  # base quality
            output_file.write(" : {} {} {} {}".format(p, b, a, q))
        output_file.write(" # {} {} NA\n".format(len(snp_list), mq))


def calculate_coverage(read_snp_map):
    read_number = len(read_snp_map)
    cov = 0
    for read_name, snp_dict in read_snp_map.items():
        snp_list = snp_dict['snp_list']
        cov += len(snp_list)

    avg = (cov / read_number) * 100
    logging.info("Average Coverage for {}: {}%".format(read_number, avg))


def main():
    parser = argparse.ArgumentParser(
        description='Convert an alignment file .bam in a .wif Hapcol compatible file')
    parser.add_argument('-b', action='store', dest='fileBAM', help='Alignment file in BAM format.', required=True)
    parser.add_argument('-vf', action='store', dest='variantFile',
                        help='SNPs variant file. Retrieved with get.variants.py script', required=True)
    parser.add_argument('-o', action='store', dest='outputDir',
                        help='Output (root) directory. Default: current directory')
    parser.add_argument('-file-delimiter', action='store', dest='fileDelim',
                        help='Set the file delimiter for the variant file. Default: \\t')
    parser.add_argument('-c', action='store_true', dest='coverageFlag',
                        help="If set calculate reads coverage and return to stdout", default=False)
    parser.add_argument('--chromosome', '-chr', action='store', dest='chromosome',
                        help='Chromosome to analyze', default='chr1')
    parser.add_argument('--base-quality', '-bq', action='store', dest='baseQuality',
                        help='Set the minimum quality for base', type=int)
    parser.add_argument('--read-quality', '-rq', action='store', dest='readQuality',
                        help='Set the minimum quality for read', type=int)
    parser.add_argument('-v', help='increase output verbosity', action='count')

    args = parser.parse_args()

    if args.v == 0 or not args.v:
        log_level = logging.INFO
    else:
        log_level = logging.DEBUG

    logging.basicConfig(level=log_level,
                        format='%(levelname)-8s [%(asctime)s]  %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    logging.info('BamToWif: Program Started')

    if args.outputDir and not os.path.exists(args.outputDir):
        logging.error('Output dir not found.')
        sys.exit(1)

    if args.outputDir:
        out_dir = args.outputDir + '/'
    else:
        out_dir = os.getcwd() + '/'

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    ann_file_delimiter = '\t'
    if args.fileDelim is not None:
        ann_file_delimiter = args.fileDelim

    logging.info('#### Calculating SNPs info ####')

    chromosome = 'chr1'
    if args.chromosome:
        chromosome = args.chromosome

    snps_map_info = load_snps(args.fileBAM, args.variantFile, ann_file_delimiter, chromosome, args.baseQuality,
                              args.readQuality)

    if args.coverageFlag:
        logging.info('#### Calculating coverage ####')
        calculate_coverage(snps_map_info)

    logging.info('#### Writing Wif ####')
    write_wif(snps_map_info, out_dir)

    logging.info('BamToWif: Program Completed')


if __name__ == '__main__':
    main()

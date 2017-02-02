usage = '''

usage: python get.variants.py prefix [variants.vcf]

given a file variants.vcf (default: stdin) containing variants and a
string prefix, output to a file prefix_i.var the variants that
correspond to contig i (e.g., a chromsome) for each unique contig i of
the #CHROM column

'''

import sys

#
# write variants to a file prefix_i.var (do nothing if i is 'None')
def write_variants(prefix, i, variants) :

    if not i : return

    filename = prefix + '_' + i + '.var'
    handle = open(filename,'w')

    for pos in sorted(variants) :
        line = str(pos) + '\t' + '\t'.join(variants[pos]) + '\n'
        handle.write(line)

    handle.close()

#
# PARSER
#

prefix = None
entree = sys.stdin

a = sys.argv[1:]
assert len(a), usage
i = 0
while i < len(a) : # bash-style argparse

    if not prefix :
        prefix = a[i]
    else :
        entree = open(a[i],'r')

    i += 1 # shift once by default in any case

#
# MAIN
#

contig = None
variants = {}
for line in entree :

    if line.startswith('#') : # burn through head
        continue

    s = line.split()[:5]
    if s[0] != contig : # new contig
        write_variants(prefix, contig, variants) # write old contig
        contig = s[0] # and reset values
        variants = {}

    # for now we take only SNVs
    if len(s[3]) == len(s[4]) == 1 :
        pos = int(s[1])
        assert pos not in variants # sanity check
        variants[pos] = s[3:] # ref and alt allele (for now)

# and write final contig
write_variants(prefix, contig, variants)

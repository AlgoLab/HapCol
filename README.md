
  HapCol
==========

A fast and memory-efficient method for haplotype assembly from long gapless
reads, like those produced by SMRT sequencing technologies (PacBio RS II) and
Oxford Nanopore flow cell technologies (MinION).

HapCol implements a fixed-parameter algorithm for the k-constrained Minimum
Error Correction problem (k-cMEC), a variant of the well-known MEC problem where
the maximum number of corrections per column is bounded by an integer k.
HapCol, while is as accurate as other exact state-of-the-art combinatorial approaches,
is significantly faster and more memory-efficient than them.
Moreover, HapCol is able to process datasets composed of both long reads (over
100 000bp long) and coverages up to 25x on standard workstations/small servers,
whereas the other approaches cannot handle long reads or coverages greater than 20x.


## Citation ##

The detailed description of the algorithm, along with an experimental comparison
with other state-of-the-art haplotype assembly tools, is presented in:

Yuri Pirola*, Simone Zaccaria*, Riccardo Dondi, Gunnar W. Klau, Nadia Pisanti, and Paola Bonizzoni  
_HapCol: Accurate and Memory-efficient Haplotype Assembly from Long Reads_.  
Bioinformatics.
*_Joint first authors_

[doi:10.1093/bioinformatics/btv495](http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btv495?ijkey=2dl7qCgbFQ9eHFj)


## Compilation ##

HapCol is distributed only on source form.
It has been developed and tested on Ubuntu Linux but should work on (or should
be easily ported to) on MacOS X.

The latest stable version of HapCol can be downloaded from GitHub in either
[.tar.gz](https://github.com/algolab/hapcol/tarball/master) or in
[.zip](https://github.com/algolab/hapcol/zipball/master) format.
Previous stable releases can be downloaded from
<https://github.com/algolab/hapcol/releases>.

HapCol depends on:

- CMake (>= 2.8)
- GNU make
- Boost program_options (tested with 1.48, previous versions should work)

We suggest to build HapCol out-of-tree with the following commands:

    mkdir -p build
    cd build
    cmake ../src
    make

The resulting file `hapcol` is the standalone executable program.

## Basic usage ##

The execution of HapCol requires to specify at least two parameters:

- `--input` (or `-i`), which specifies the file containing the reads in input (in
  WIF format);
- `--output` (or `-o`), which specifies the file for the computed haplotypes.

Optional parameters are:

- `--error-rate` (or `-e`), for specifying the estimated sequencing error rate
  of the input reads;
- `--alpha` (or `-a`), for specifying the significance level (lower levels
  require more computational resources but increase the probability of finding a
  feasible solution);
- `--discard-weights` (or `-u`), for discarding weights while computing the
  optimal solution (notice that accuracy of the reconstructed haplotypes may
  decrease);
- `--no-ambiguous` (or `-x`), do not mask ambiguous positions in the output
  haplotypes with a `X`.
- `--unique` (or `-U`), do not split the input into independent blocks and
  consider the input as a unique block (notice that by default the input is
  split into independent blocks and, in correspondence to these, brackets
  are added to the output).
- `--all-heterozygous` (or `-A`), for solving the input instance under the
  traditional all-heterozygous assumption.

For example, HapCol can be executed on the sample data included with the program
with the following command (given from the directory `build/`):

    ./hapcol -i ../docs/sample.wif -o haplotypes.txt

which should save a solution of cost 62 in the weighted case (or cost 7 in
the unweighted case, if flag `-u` is added) in file `haplotypes.txt`.


## Data ##

The repository contains in the folde `data` some of the main simulated datasets and the real dataset that have been used in (Pirola, Zaccaria et al., Bioinformatics, 2016).

### Input format ###

The current version of HapCol assumes that we already selected interesting positions (such as SNPs) and takes as input a dataset in WIF format.
The WIF format (introduced in (Patterson et al., RECOMB, 2014)) is defined as follows:

- Each line corresponds to a read and read must be ordered depedning on the starting position. Each line terminated with a string preeceded by the comment symbol `#` (such as `: # 0 0 : N N`) that is not currently considered
- Each read is composed of the corresponding SNP alleles that are reported using `:` as separator.
- Each SNP allele is reported in the format `POS ALL BIN WEI` such that 
---- `POS` is the position of the SNP corresponding to its genomic coordinate (or column id in the fragment matrix)
---- `ALL` is the allele of the corresponding SNP, equal to a value in {A, C, G, T}
---- `BIN` is the binary representation of the alle such that 0 correspons to the MAJOR ALLELE (reference allele) and 1 corresponds to the MINOR ALLELE (variant)
---- `WEI` is the phred score of the SNP reported in the numeric format

Here in the following an example of a wif format:

    1 A 0 30 :  # 0 0 : N N
    1 G 1 33 : 2 A 0 31 : # 0 0 : N N
    1 G 1 29 : 3 T 1 30 : # 0 0 : N N
    2 A 0 31 : # 0 0 : N N

This example is composed of 4 reads. The first read covers only position 1 having a major allele with a phred score of 30.
The second reads covers positions 1 and 2 where the first has a minor allele with score of 33 and the second has a major allele with score of 31.
The third reads covers positions 1 and 3 having a minor allele and score of 29 in the first and a minor allele and score of 30 in the third.
Note that the third read does not cover position 2 and the corresponding gap is not indicated in this format.
Last, the fourth reads covers only position 1 with a major allele and score of 31. 


### Simulated ###

The simulated datasets are located in `data/simulated`. Each simulated dataset has a file name that encodes the main characteristics of its data:

     Venter-chr$C.artificial$A.length$L.cov$V.single.shuffled$S.H$H.00.wif

where the parameters in the name have the following meaning:

PARAMETER|
---------|-------------------------------
$C       | The corresponding chromosome that is either 1 or 15
$A       | The substitution-error rate that is either 5 (%) or 1 (%)
$L       | The length of the original simulated reads (5000, 10000, 50000)
$V       | The original coverage of the simulated reads that is always 30 in these datasets
$S       | The random seeds used to randomly select a subset of the reads by downsampling
$H       | The maximum coverage of the current dataset

Moreover the datasets are subdivided into two folders:

- `10-indel-errorrate` contains the datasets with a indel-error rate of 10%

- `15-indel-errorrate` contains the datasets with a indel-error rate of 15%

Last, folder `true-haplotypes` contains the true haplotypes used for the simulations. In particular, for each chromosome the haplotypes are given in a string format: each file (`Venter-chr1.haplo` and `Venter-chr15.haplo`) contains two lines, each corresponding to a haplotype of the related chromosome.
For each chromosome the coresponding file `Venter-chr1.positions` and `Venter-chr15.positions`, respectively, contains the genomic coordinates of each SNP in the haplotypes.
In the latter, each line of the file corresponds to a different SNP and they are reported in the same order as they appear in the haplotype strings.

### Real ###

The real dataset is the one published in (Duitama et al., NAR, 2012) released at https://owww.molgen.mpg.de/~genetic-variation/SIH/data/.
The folder `input-duitama` contains all the datasets, one for each chromosome.
In particular, the subfolder `original` contains the original datasets, instead `wif` contains the corresponding dataset of each chromosome in wif format (ready to execute with hapcol).
Last, `haplotypes-duitama` contains the phase of the corresponding SNPs: For each chromosome, the corresponding file contains a line for each SNP such that the first column indicates its genomic coordinate, the second indicates the allele of the first haplotype, and, hence, the last column indicates the allele on the second haplotype.



## License ##

HapCol is licensed under the terms of the GNU GPL v2.0


## Contacts ##

For questions or support, please contact <simone.zaccaria@disco.unimib.it>
or <yuri.pirola@disco.unimib.it>

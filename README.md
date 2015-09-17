
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

Yuri Pirola, Simone Zaccaria, Riccardo Dondi, Gunnar W. Klau, Nadia Pisanti, and Paola Bonizzoni  
_HapCol: Accurate and Memory-efficient Haplotype Assembly from Long Reads_.  
Bioinformatics.
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



## License ##

HapCol is licensed under the terms of the GNU GPL v2.0


## Contacts ##

For questions or support, please contact <simone.zaccaria@disco.unimib.it>
or <yuri.pirola@disco.unimib.it>

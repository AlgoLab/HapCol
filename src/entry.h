/**
 *
 * This file was originally part of
 * WhatsHap (https://bitbucket.org/whatshap/whatshap), written by:
 * Tobias Marschall <t.marschall@mpi-inf.mpg.de>,
 * Marcel Martin <marcel.martin@scilifelab.se>,
 * Murray Patterson,
 * Alexander Schönhuth
 *
 * and released under the following licence:
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * The file was partially modified by Yuri Pirola and Simone Zaccaria
 *
 **/

#ifndef ENTRY_H
#define ENTRY_H

#include <iostream>

typedef long int readid_t;

#define SENTINEL_READID (-1)

class Entry {

public:

  typedef enum { MAJOR_ALLELE = 0, MINOR_ALLELE = 1, BLANK = 2, EQUAL_SCORES = 3 } allele_t;

  Entry(readid_t r, allele_t m, unsigned int p)
      :read_id(r), allele_type(m), phred_score(p)
  {}

  readid_t get_read_id() const { return read_id; }
  allele_t get_allele_type() const { return allele_type; }
  unsigned int get_phred_score() const { return phred_score; }

  void set_read_id(readid_t r) { read_id = r; }
  void set_allele_type(allele_t m) { allele_type = m;}
  void set_phred_score(unsigned int p) { phred_score = p; }

  friend std::ostream& operator<<(std::ostream& out, const Entry& e);

private:
  readid_t read_id;
  allele_t allele_type;
  unsigned int phred_score;
};

#endif

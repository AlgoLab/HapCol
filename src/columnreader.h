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
 **/

/*
  reads and outputs, incrementally, the columns of a .wif file.
  ColumnReader will never output a column of height larger than
  'coverage_threshold', even if it reads in a column of height larger
  than this (in which case discards columns until height is within
  this threshold, according to some critieria) -- Murray Patterson,
  Sept 2013

  obacht : think of some sort of criteria.  For now, threshold is a
  cutoff (that is, it assumes that the preprocessing step took care of
  controlling the coverage) -- murray
*/

#ifndef COLUMN_READER_H
#define COLUMN_READER_H

#include "entry.h"
#include <string>
#include <fstream>
#include <vector>
#include <queue>
#include <list>
#include <memory>

class ColumnReader {

public:
  /** Constructor.
   *  @param remove_weights If true, all weights (i.e. phred scores) are removed,
   *                        that is, set to one. */
  ColumnReader(const std::string& f, const unsigned int c, const bool remove_weights = false);

  // destructor
  ~ColumnReader() {
    // Nothing to do
  }

  // number of columns (snp positions)
  size_t num_cols() { return positions.size(); }

  // number of rows (so far)
  size_t num_rows() { return row; }


  // true iff the reader has another column
  bool has_next();
  // get next column
  std::vector<Entry> get_next();
  // return a const pointer to the positions
  const std::vector<unsigned int>& get_positions() const {
    return positions;
  }

private:
  std::ifstream ifs;
  const unsigned int coverage_threshold;
  const bool remove_weights_;

  std::vector<unsigned int> positions; // entry (snp) positions (indexed by column)
  typedef std::queue <Entry> queue_t;
  typedef std::list<queue_t> buffer_t;
  buffer_t buffer; // buffer on file

  size_t row; // current row (column) during the process
  size_t column;

  // auxiliary functions

  /*
    reads in .wif file, gets its entry (snp) positions and places
    them positions vector and returns true iff all of this was
    successful
  */
  bool compute_positions();

  /*
    returns true if buffer is not leftmost total (a buffer is
    'leftmost total' when all start positions of rows (row suffixes)
    in the buffer are the same (as current column), and does not
    read in this case, o.w. it reads line from .wif file into buffer
    and returns true if this causes the buffer to no longer be
    leftmost total, i.e., when the start position of this line is
    beyond the current column)

    precondition: the file stream has 'good' status
  */
  bool read_line();

};

std::ostream& operator<<(std::ostream& out, const std::vector<unsigned int>& v);

#endif

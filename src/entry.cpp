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

#include <cassert>

#include "entry.h"

std::ostream& operator<<(std::ostream& out, const Entry& e) {
  out << "Entry(" << e.read_id << ',';
  switch (e.allele_type) {
    case Entry::MAJOR_ALLELE:
      out << "MAJOR";
      break;
    case Entry::MINOR_ALLELE:
      out << "MINOR";
      break;
    case Entry::BLANK:
      out << "BLANK";
      break;
    case Entry::EQUAL_SCORES:
      out << "EQUAL_SCORES";
      break;
    default:
      assert(false);
  }
  out << ',' << ((int)e.phred_score) << ')';
  return out;
}

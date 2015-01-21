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

#include "entry.h"
#include "columnreader.h"
#include <algorithm> // for 'sort'
#include <exception>
#include <stdexcept>
#include <set>
#include <cassert>

using namespace std;

// constructor
ColumnReader::ColumnReader(const string& f, const unsigned int c, const bool remove_weights)
  :   ifs(f.c_str()), coverage_threshold(c), remove_weights_(remove_weights) {

  const bool s = compute_positions(); // get the positions
  if (!s) throw runtime_error("error when reading positions");
  ifs.close();

  ifs.open(f.c_str(), ios::in); // open again for the rest
  row = 0;
  column = 0;
}

// true iff reader (buffer) has a next column
bool ColumnReader::has_next() {

  // keep reading lines until start position passes current column
  while(ifs.good()) { if(read_line()) { break; }}
  
  if(buffer.empty()) { return false; }
  return true;
}

// outputs next column from buffer, thresholding if necessary
vector<Entry> ColumnReader::get_next() {
  assert(!buffer.empty());
//     cerr << "ColumnReader::get_next(): start" << endl;
  
  vector<Entry> c;

  unsigned int t = 0; // for now, threshold is just a cutoff
  buffer_t::iterator it = buffer.begin();
  while (it != buffer.end()) {
    queue_t& q = *it;
    assert(!q.empty());
//       cerr << "  considering next queue, first entry: ";
//       if (q.front()==0) cerr << "DUMMY" << endl;
//       else cerr << (*q.front()) << endl;
    if (t > coverage_threshold) {
      throw runtime_error("Error: coverage threshold exceeded!");
    }
    if(q.front().get_read_id() != -1) { // entry not a dummy
//         cerr << "    --> pushing entry" << endl;
      c.push_back(q.front());
      ++t;
    }
    q.pop(); // pop queue, and if empty now then erase it
    if(q.empty()) {
//         cerr << "    queue now empty, removing it" << endl;
      it = buffer.erase(it);
    } else {
      ++it;
    }
  }
  column += 1;
//     cerr << "ColumnReader::get_next(): end" << endl;
  return c;
}

// auxiliary functions

/*
  reads in .wif file, gets its entry (snp) positions and places them
  positions vector and returns true iff all of this was successful
*/
bool ColumnReader::compute_positions() {

  set<unsigned int> ps; // for adding, incrementally, (many duplicate) positions
  if(ifs.good()) {

    unsigned int p; // each line (read) has at least one entry
    ifs >> p;
    bool eol = false; // 'end of line'
    while(ifs.good()) {
      eol = false; // for 2 or more lines
      ps.insert(p); // insert (already read) position

      string s; // temporary input holder
      ifs >> s; ifs >> s; ifs >> s; // discard everything else in entry

      // next entry ...
      ifs >> s; // discard ":"
      ifs >> s; // either a position, "--" (gap) or "#" (end of line)
      if(s == "--") { ifs >> s; ifs >> s; } // discard ":", store next

      if(s == "#") { eol = true;
        ifs >> s; // fragment map quality
        ifs >> s; // either another (paired end) mapq, or ":"
        if(s == ":") { ifs >> s; } // discard ":" and 'u' character
        else { ifs >> s; ifs >> s; ifs >> s; } // discard ":" and both 'u' chars
        if(ifs.good()) { ifs >> p; }
      } // get position of next line
      else { 
        p = atoi(s.c_str());
      } // else it was another entry
    } 

    if(!eol) { // file ended before end of line
      throw runtime_error("error : line has no terminator");
      return false; 
    }
  } else { 
    return false;
  } // something was wrong with the file

  // now, feed set into positions vector and sort
  
  set<unsigned int>::iterator it;
  for(it=ps.begin(); it!=ps.end(); ++it) {
    positions.push_back(*it);
  }

  sort(positions.begin(), positions.end());

  return true;
}

/*
  returns true if buffer is not leftmost total (a buffer is
  'leftmost total' when all start positions of rows (row suffixes)
  in the buffer are the same (as current column), and does not read
  in this case, o.w. it reads line from .wif file into buffer and
  returns true if this causes the buffer to no longer be leftmost
  total, i.e., when the start position of this line is beyond the
  current column)

  precondition: the file stream has 'good' status
*/
bool ColumnReader::read_line() {
//     cerr << "ColumnReader::read_line(): start" << endl;
  // if buffer is not leftmost total, return true, else continue
  if(!buffer.empty()) { // (an empty buffer is leftmost total, by definition)
    if((*buffer.rbegin()).front().get_read_id() == -1) { 
//         cerr << "ColumnReader::read_line(): end (leftmost total already)" << endl;
      return true; 
    }
  }
  
  bool a = false; // start position is larger than current column
  size_t c = column; // current column as we read this line
  queue<Entry> q; // queue that stores the line

  unsigned int p; // each line (read) has at least one entry
  ifs >> p;
  // is the file over?
  if (!ifs.good()) {
    return true;
  }
  while(p > positions[c]) {
//       cout << "  Pushing 0 entry" << endl;
    q.push(Entry(-1, Entry::BLANK, -1)); // push dummy with left padding (-2) tag
    ++c;
    a = true;
  }
  
  bool eol = false; // 'end of line'
  while(ifs.good()) {

    if (eol) { break; }

    while(p > positions[c]) { // to deal with gaps
//         cout << "  Pushing (-1)" << endl;
      q.push(Entry(row, Entry::BLANK, 0)); // fill gap
      ++c;
    }

    bool m;
    unsigned int ph;
    string s; // temporary input holder
    
    ifs >> s; // discard the nucleotide
    ifs >> m;
    ifs >> ph;
//       cout << row << " " << m << " " << ph << " : ";
//       cout << "  Pushing " << Entry(row, (m?MINOR_ALLELE:MAJOR_ALLELE), ph) << endl;
    if (remove_weights_) {
      q.push(Entry(row, (m?Entry::MINOR_ALLELE:Entry::MAJOR_ALLELE), 1));
    } else {
      q.push(Entry(row, (m?Entry::MINOR_ALLELE:Entry::MAJOR_ALLELE), ph));
    }
    ++c;
    
    // next entry ...
    ifs >> s; // discard ":"
    ifs >> s; // either a position, "--" (gap) or "#" (end of line)
    if(s == "--") { ifs >> s; ifs >> s; } // discard ":", store next
    
    if(s == "#") {
      eol = true;
      ifs >> s; // fragment map quality
      ifs >> s; // either another (paired end) mapq, or ":"
      if(s == ":") { ifs >> s; } // discard ":" and 'u' character (for now)
      else { ifs >> s; ifs >> s; ifs >> s; }
    } // discard ":" and both 'u' chars
    else { 
      p = atoi(s.c_str());
    }
  }

  if(!eol) { // file ended before end of line
    throw runtime_error("error : line has no terminator");
  }

  // add queue, update row and return
  buffer.push_back(q);
  ++row;
//     cout << endl << "a is : " << a << endl;
//     cerr << "ColumnReader::read_line(): end" << endl;
  return a;
}

std::ostream& operator<<(std::ostream& out, const std::vector<unsigned int>& v) {
  std::vector<unsigned int>::const_iterator it = v.begin();
  for (; it != v.end(); ++it) {
    if (it != v.begin()) out << ',';
    out << (*it);
  }
  return out;
}
